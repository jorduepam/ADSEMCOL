using Test

include("MemoryReporting.jl")
using .MemoryReporting

@testset "Memory reporting" begin
    withenv("ADSEMCOL_MEMORY_REPORT_STEPS" => "0") do
        @test MemoryReporting.reporter() === nothing
    end
    for value in ("-1", "invalid")
        withenv("ADSEMCOL_MEMORY_REPORT_STEPS" => value) do
            @test_throws ErrorException MemoryReporting.reporter()
        end
        withenv("SLURM_JOB_ID" => nothing) do
            @test MemoryReporting.job_cgroup() === nothing
        end
    end

    mktempdir() do job
        step = joinpath(job, "step_0")
        task = joinpath(step, "task_0")
        mkpath(task)
        gib = UInt64(1024)^3
        write(joinpath(job, "memory.current"), string(10gib))
        write(joinpath(job, "memory.peak"), string(15gib))
        write(joinpath(job, "memory.max"), string(40gib))
        write(joinpath(step, "memory.current"), string(9gib))
        write(joinpath(step, "memory.max"), string(32gib))
        write(joinpath(task, "memory.max"), "max")

        snapshot = MemoryReporting.cgroup_memory(job, task)
        @test snapshot.current == 10gib
        @test snapshot.peak == 15gib
        @test snapshot.limit == 32gib
        @test snapshot.limit_cgroup == step
        @test snapshot.limited_current == 9gib

        reporter = MemoryReporting.Reporter(500, job, task)
        mktemp() do _, output
            redirect_stdout(output) do
                MemoryReporting.report(reporter, 499)
                MemoryReporting.report(reporter, 500)
            end
            seekstart(output)
            text = read(output, String)
            @test length(split(chomp(text), '\n')) == 1
            @test occursin("Memory step=500 Julia_peak_RSS=", text)
            @test occursin("job_current=10.000 GiB job_peak=15.000 GiB", text)
            @test occursin("limit_scope=step_0 scope_current=9.000 GiB limit=32.000 GiB", text)
        end

        write(joinpath(job, "memory.max"), "max")
        write(joinpath(step, "memory.max"), "max")
        rm(joinpath(job, "memory.peak"))
        snapshot = MemoryReporting.cgroup_memory(job, task)
        @test snapshot.limit === nothing
        @test snapshot.peak === nothing

        write(joinpath(job, "memory.current"), string(41gib))
        write(joinpath(job, "memory.peak"), string(42gib))
        mktemp() do _, output
            redirect_stdout(output) do
                MemoryReporting.report(reporter, 500)
            end
            seekstart(output)
            text = read(output, String)
            @test occursin("job_current=41.000 GiB job_peak=42.000 GiB", text)
            @test occursin("job_step_task_limit=unlimited", text)
        end

        rm(joinpath(step, "memory.max"))
        mktemp() do _, output
            redirect_stdout(output) do
                MemoryReporting.report(reporter, 500)
            end
            seekstart(output)
            @test occursin("job_cgroup=unavailable", read(output, String))
        end
    end
end
