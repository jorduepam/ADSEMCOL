if get(ENV, "ADSEMCOL_REPORT_MAXRSS", "0") == "1"
    atexit(() -> println("Peak Julia process RSS: $(Sys.maxrss()) bytes"))
end

#Load problem:
ProblemName     = ARGS[1]
include("$(ProblemName).jl")

try
    if ProblemName=="DetonationWave"
        tf              = parse(Float64, ARGS[2])
        delta           = parse(Float64, ARGS[3])
        epsilon         = parse(Float64, ARGS[4])
        TolS            = parse(Float64, ARGS[5])
        TolT            = parse(Float64, ARGS[6])
        Deltat_SaveRes  = parse(Float64, ARGS[7])
        SC              = parse(Int, ARGS[8])
        DetonationWave(tf=tf, delta=delta, epsilon=epsilon, 
            TolS=TolS, TolT=TolT, 
            SaveRes=true, Deltat_SaveRes=Deltat_SaveRes, SC=SC)
            
    elseif ProblemName=="ShockFlame_FXP"
        FesOrder        = parse(Int, ARGS[2])
        tf              = parse(Float64, ARGS[3])
        q               = parse(Float64, ARGS[4])
        beta            = parse(Float64, ARGS[5])
        ML              = parse(Float64, ARGS[6])
        epsilon         = parse(Float64, ARGS[7])
        Lx              = parse(Float64, ARGS[8])
        SC              = parse(Int, ARGS[9])
        nu              = parse(Float64, ARGS[10])
        kappa_rho_cv    = parse(Float64, ARGS[11])
        ShockFlame_FXP(FesOrder=FesOrder, tf=tf, q=q, 
            RThatI=beta, RThatB=beta, ML=ML, epsilon=epsilon, 
            Lx=Lx, kappa_rho_cv=kappa_rho_cv, nu=nu,
            SaveRes=true, Deltat_SaveRes=tf/50, SC=SC)

    elseif ProblemName=="ShockFlame2_FXP"
        FesOrder        = parse(Int, ARGS[2])
        tf              = parse(Float64, ARGS[3])
        Lx              = parse(Float64, ARGS[4])
        q               = parse(Float64, ARGS[5])
        beta            = parse(Float64, ARGS[6])
        ML              = parse(Float64, ARGS[7])
        epsilon         = parse(Float64, ARGS[8])
        nu              = parse(Float64, ARGS[9])
        kappa_rho_cv    = parse(Float64, ARGS[10])
        D               = parse(Float64, ARGS[11])
        SC              = parse(Int, ARGS[12])
        ShockFlame2_FXP(FesOrder=FesOrder, tf_hat=tf, Lx_hat=Lx, q=q,
            RThatI=beta, RThatB=beta, ML=ML, epsilon_hat=epsilon,
            nu_hat=nu, kappa_rho_cv_hat=kappa_rho_cv, D_hat=D,
            SaveRes=true, Deltat_SaveRes=tf/100/50, SC=SC)
            
    else
        error("Undefined problem $(ProblemName)")
    end
    
catch err
    
    println(Error2String(err))
    rethrow()
            
end
