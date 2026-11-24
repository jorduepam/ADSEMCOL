module MemoryReporting

using Printf

struct Reporter
    interval::Int
    job_cgroup::Union{Nothing,String}
    task_cgroup::Union{Nothing,String}
end

function job_cgroup()
    job_id = get(ENV, "SLURM_JOB_ID", "")
    isempty(job_id) && return nothing
    for line in eachline("/proc/self/cgroup")
        startswith(line, "0::/") || continue
        parts = split(line[5:end], '/'; keepempty=false)
        job_index = findlast(==("job_$job_id"), parts)
        job_index === nothing && return nothing
        job = joinpath("/sys/fs/cgroup", parts[1:job_index]...)
        task = joinpath("/sys/fs/cgroup", parts...)
        isdir(job) && isdir(task) && return (job, task)
    end
    return nothing
end

function reporter()
    value = get(ENV, "ADSEMCOL_MEMORY_REPORT_STEPS", "0")
    interval = tryparse(Int, value)
    (interval !== nothing && interval >= 0) ||
        error("ADSEMCOL_MEMORY_REPORT_STEPS must be a nonnegative integer: $value")
    interval == 0 && return nothing

    cgroup = try
        job_cgroup()
    catch err
        (err isa Base.IOError || err isa SystemError) || rethrow()
        println("Job cgroup memory unavailable: $(sprint(showerror, err))")
        nothing
    end
    if cgroup === nothing
        println("Job cgroup memory unavailable: no accessible Slurm cgroup v2 found.")
        flush(stdout)
        return Reporter(interval, nothing, nothing)
    end
    return Reporter(interval, cgroup...)
end

function cgroup_memory(job::String, task::String)
    current = parse(UInt64, strip(read(joinpath(job, "memory.current"), String)))
    peak_file = joinpath(job, "memory.peak")
    peak = isfile(peak_file) ? parse(UInt64, strip(read(peak_file, String))) : nothing

    # Check job, step and task limits, but not shared ancestors outside the job.
    limit = nothing
    limit_cgroup = job
    path = task
    while true
        raw = strip(read(joinpath(path, "memory.max"), String))
        value = raw == "max" ? nothing : parse(UInt64, raw)
        if value !== nothing && (limit === nothing || value < limit)
            limit = value
            limit_cgroup = path
        end
        path == job && break
        parent = dirname(path)
        parent == path && error("Task cgroup is not within job cgroup")
        path = parent
    end
    limited_current = limit === nothing ? nothing :
        parse(UInt64, strip(read(joinpath(limit_cgroup, "memory.current"), String)))
    return (; current, peak, limit, limit_cgroup, limited_current)
end

function report(reporter::Reporter, step::Int)
    step > 0 && step % reporter.interval == 0 || return
    gib = 1024.0^3
    @printf("Memory step=%d Julia_peak_RSS=%.3f GiB", step, Sys.maxrss() / gib)
    if reporter.job_cgroup !== nothing
        try
            memory = cgroup_memory(reporter.job_cgroup, reporter.task_cgroup)
            @printf(" job_current=%.3f GiB", memory.current / gib)
            if memory.peak === nothing
                print(" job_peak=unavailable")
            else
                @printf(" job_peak=%.3f GiB", memory.peak / gib)
            end
            if memory.limit === nothing
                print(" job_step_task_limit=unlimited")
            else
                scope = relpath(memory.limit_cgroup, reporter.job_cgroup)
                @printf(" limit_scope=%s scope_current=%.3f GiB limit=%.3f GiB",
                    scope == "." ? "job" : scope, memory.limited_current / gib, memory.limit / gib)
            end
        catch err
            (err isa Base.IOError || err isa SystemError || err isa ArgumentError) || rethrow()
            print(" job_cgroup=unavailable ($(sprint(showerror, err)))")
        end
    else
        print(" job_cgroup=unavailable")
    end
    println()
    flush(stdout)
end

end
