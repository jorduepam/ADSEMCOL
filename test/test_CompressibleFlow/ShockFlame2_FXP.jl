include("test_CompressibleFlow.jl")
include("InputData/ReactiveGasAuxFunctions.jl")

"""
    ShockFlame2_FXP(; kwargs...)

Simulate the interaction of a flame and a shock with the reactive `GasFXP`
model on an adaptive two-dimensional mesh. The initial flame and shock are
at `10*delta_flame` and `20*delta_flame`; the domain extends from zero to
`Lx` and has height `delta_flame`. The simulation initializes the solution,
advances it in time, optionally plots and saves results, and returns the solver.
It stops early if a time step reports a nonpositive convergence flag.

# Model and scaling

The fresh gas has `p0 = RT0 = 1`, `rho0 = p0/RT0` and `gamma = 1.4`.
The reference scales are `uc = sqrt(RT0)`, `tc = 1/BR`, and `Lc = uc*tc`.
Parameters ending in `_hat` are dimensionless and converted to solver units:
`tf_hat`, `Lx_hat`, `delta_flame_hat`, and `delta_shock_hat` are scaled by
`tc` or `Lc`, while `epsilon_hat`, `nu_hat`, `kappa_rho_cv_hat`, and `D_hat`
are multiplied by `Lc^2/tc`.

- `q`, `ML`: heat-release parameter and Mach number used for the initial
  states. `YF0`, `YX0`: initial fresh-gas fuel and intermediate mass fractions.
- `RThatI`, `RThatB`: activation-temperature ratios relative to the
  Chapman--Jouguet temperature. `BIhat`, `BBhat`, `BR`: reaction-rate
  parameters used to construct the `GasFXP` kinetic coefficients.
- `SpoungeBool` enables the sponge source; `SpoungeLc_hat` is its length
  scale as a fraction of `Lx`. The commented block after `SpoungeFun` plots
  the actual damping coefficient versus `x/Lx`, even when the source is
  disabled. The Gaussian has no sharp start at `Lx - SpoungeLc_hat*Lx`.

# Numerical and output options

- `FesOrder`, `RKMethod`: finite-element order and time-integration method.
  `KC35` selects `IRK_Step!`; other methods use `LIRKHyp_Step!`.
- `TolS`, `AMA_MaxIter`, `AMA_SizeOrder`, `AMA_AnisoOrder`: spatial
  adaptation tolerance and mesh-adaptation controls. `TolT` and `Deltat0`
  set the time tolerance and initial time step; `Deltat0` is in solver time.
- `PlotFig` enables plotting; `PlotVars`, `MarkSF`, `wFig`, `hFig`, `mFig`
  and `nFig` control the displayed variables, markers and figure layout.
  `Nt_SaveFig` and `Deltat_SaveFig` control the plotting interval in steps
  and solver time. `SaveFig` writes figures only when plotting is enabled.
- `SaveRes` writes numbered JLD2 snapshots for study case `SC` at the
  intervals `Nt_SaveRes` and `Deltat_SaveRes` (steps and solver time). Each
  snapshot includes the solver, model parameters and histories `tv`, `xsv`,
  `xfv` (in solver units), and dimensionless `Mrelv`. The final `_info.jld2`
  is written when the final time is reached. Output locations come from
  `ResUbi.jl`.
- `CSS`, `cmap`, and `LFact` are accepted by this function but are not used
  by its current implementation.

# Front positions

At initialization and after each accepted step, pressure and heat release
`Q_reac` are interpolated along the domain centre line at the distinct mesh
x coordinates. The flame position is the sample with maximum `Q_reac` on
that line, without sub-cell refinement; if no sample is finite it is `NaN`.

The shock position is the first positive local maximum of the pressure
gradient scanning from fresh gas at smaller x, with pressure at the peak
sample strictly greater than `(1 + shock_alpha)*p0`. `shock_alpha` is
dimensionless, positive, and defaults to `0.3`; for example, `0.1` gives
`1.1*p0`. The gradient accounts for nonuniform spacing. Isolated peaks are
refined with a three-point parabola; a flat peak uses its middle sample.
No causal filter is applied. The threshold is not a minimum total pressure
jump: weak shocks can be missed and other waves or noise can pass it.

If no shock peak qualifies, a warning is emitted and its position is `NaN`.
The shock speed and relative Mach number derived from it are also `NaN` for
that step and the following step. Saved snapshots include `shock_alpha` and
`p0`. The launcher uses the default unless the function is called directly.
"""
function ShockFlame2_FXP(; FesOrder::Int=3, 
    tf_hat::Float64=20.0, RKMethod::String="Ascher3",
    epsilon_hat::Float64=0e-2, nu_hat::Float64=1e-2, kappa_rho_cv_hat::Float64=0e-6, D_hat::Float64=1e-2, 
    YF0::Float64=1.0-0e-8, YX0::Float64=0.0,
    delta_flame_hat::Float64=1.0, delta_shock_hat::Float64=max(epsilon_hat,nu_hat),
    q::Float64=5.0, ML::Float64=2.0, shock_alpha::Float64=0.3,
    RThatI::Float64=10.0, RThatB::Float64=10.0, #RThat=RT_a/RT_1, #1: means post-shock condition
    BIhat::Float64=2.5e-3, BBhat::Float64=4.4e5, BR::Float64=100.0, 
    Lx_hat::Float64=300.0, SpoungeBool::Bool=true, SpoungeLc_hat::Float64=0.1, 
    #
    TolS::Float64=1e-6, AMA_MaxIter::Int=200, AMA_SizeOrder::Int=FesOrder, AMA_AnisoOrder::Int=2, 
    #
    TolT::Float64=1e-4, Deltat0::Float64=1e-6,
    #
    PlotFig::Bool=false, MarkSF::Bool=true, wFig::Float64=9.50, hFig::Float64=6.50, 
    PlotVars::Vector{String}=String[], 
    SaveFig::Bool=false, Nt_SaveFig::Int=5, Deltat_SaveFig::Float64=Inf,
    mFig::Int=max(1,length(PlotVars)), nFig::Int=Int(ceil(length(PlotVars)/mFig)), cmap::String="jet",
    #
    SaveRes::Bool=false, Nt_SaveRes::Int=typemax(Int), Deltat_SaveRes::Float64=0.1, 
    #
    CSS::Float64=0.1, 
    #
    SC::Int=0, 
    LFact::Float64=1.0)     #LFact is for debugging

    isfinite(shock_alpha) && shock_alpha > 0 ||
        throw(ArgumentError("shock_alpha must be finite and positive"))

    #---------------------------------------------------------------------
    #PROBLEM DATA:
    
    #Fuel state:
    gamma           = 1.4
    p0              = 1.0
    RT0             = 1.0
    rho0            = p0/RT0
    Q0              = q*gamma/(gamma-1)*RT0/YF0
    
    #Burned gas state:
    RTL             = (1.0+q)*RT0
    pL              = p0
    rhoL            = pL/RTL
    YFL             = 0.0
    
    #Conditions at the right:
    uL              = ML*sqrt(gamma*RTL)
    rhoR            = rhoL / ( (gamma-1.0)/(gamma+1.0) + 2.0/(gamma+1.0)/ML^2 )
    uR              = rhoL*uL/rhoR
    pR              = pL+rhoL*uL^2-rhoR*uR^2
    YFR             = 0.0
    u0              = uL
    RTR             = pR/rhoR
    
    #Chapman--Jouguet state with same incoming conditions:
    Min, rhoin, uin, pin, rhoCJ, uCJ, pCJ, _ = ChapmanJouguet(gamma,q)
    rhoCJ           *= rho0/rhoin
    pCJ             *= p0/pin
    RTCJ            = pCJ/rhoCJ
    uCJ             *= sqrt(p0/rho0) / sqrt(pin/rhoin)
    
    #Characteristic velocity, time and length:
    uc              = sqrt(RT0)
    tc              = 1.0/BR
    Lc              = uc*tc
    
    tf = tf_hat * tc
    
    # Plotting scale factors (convert from physical to adimensional units for visualization)
    scale_x_plot    = 1.0 / Lc  # spatial scale: physical length → adimensional
    scale_t_plot    = 1.0 / tc  # temporal scale: physical time → adimensional


    #Define gas:
    GasModel                = GasFXP()
    GasModel.epsilon        = epsilon_hat * Lc^2 / tc
    GasModel.gamma          = gamma
    GasModel.nu             = nu_hat * Lc^2 / tc
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = kappa_rho_cv_hat * Lc^2 / tc
    GasModel.hfF            = 0.0
    GasModel.hfX            = 0.0
    GasModel.hfP            = -Q0
    GasModel.D              = D_hat * Lc^2 / tc
    GasModel.BR             = BR
    GasModel.BI             = BIhat * GasModel.BR 
    GasModel.BB             = BBhat * GasModel.BR/rhoCJ
    GasModel.RTI            = RThatI * RTCJ
    GasModel.RTB            = RThatB * RTCJ
    
    #Define domain:
    delta_flame         = delta_flame_hat * Lc
    delta_shock         = delta_shock_hat * Lc
    Lx                  = Lx_hat * Lc
    x2                  = Lx
    y2                  = 1.0*delta_flame
#     y2                  = 0.1*Lx
    x0_flame            = 10.0*delta_flame
#     x0_flame            = -10.0*Lx
    x0_shock            = 20.0*delta_flame
    MeshFile            = "$(@__DIR__)/../../temp/SC$(SC).geo"
    TrMesh_Rectangle_Create!(MeshFile, 0.0, x2, 5, 0.0, y2, 2)

    
    #Initial conditions:
    function u0fun(x::Vector{Matrix{Float64}})
        
        rhoYF           = @mlv SmoothHeaviside(x[1]-x0_flame, delta_flame, rho0*YF0, 0.0)
        rhoYX           = @mlv SmoothHeaviside(x[1]-x0_flame, delta_flame, rho0*YX0, 0.0)
        rhoYP           = @mlv SmoothHeaviside(x[1]-x0_flame, delta_flame, rho0*(1-YF0-YX0), rhoL) + 
                                SmoothHeaviside(x[1]-x0_shock, delta_shock, 0.0, rhoR-rhoL)
        rhovx           = @mlv SmoothHeaviside(x[1]-x0_flame, delta_flame, rho0*u0, rhoL*uL) + 
                                SmoothHeaviside(x[1]-x0_shock, delta_shock, 0.0, rhoR*uR-rhoL*uL)
        rhovy           = @mlv 0.0*x[1]
        rhoE            = @mlv (SmoothHeaviside(x[1]-x0_flame, delta_flame, p0, pL) + 
                                SmoothHeaviside(x[1]-x0_shock, delta_shock, 0.0, pR-pL)) / (gamma-1.0) + 
                                0.5*rhovx^2/(rhoYF+rhoYX+rhoYP) - rhoYP*Q0
        return [rhoYF, rhoYX, rhoYP, rhovx, rhovy, rhoE]
        
    end
    
    
    #Boundary conditions:
    BC_walls        = SlipAdiabatic()
    function pOutlet(t::Float64, x::Vector{Matrix{Float64}})
        
        p           = @mlv 0.0*x[1]+pR
        return [p]
        
    end
    function uInlet(t::Float64, x::Vector{Matrix{Float64}})
        
        return u0fun(x)
        
    end
    BC_outlet       = SubsonicOutlet1(FWt11((t,x)->pOutlet(t,x)))
#     BC_outlet       = SupersonicOutlet1()
#     BC_outlet       = DoNothing1()
    BC_inlet        = SupersonicInlet1(FWt11((t,x)->uInlet(t,x)))
    
    #Spounge boundary condition:
    #Add a source term Q_alpha = - sigma(t,x) * (u_alpha - uBC(t))
    #Parameters (all DIMENSIONAL):
    # - SpoungeLc_hat: Adimensional length scale of spounge zone (fraction of Lx)
    # - SpoungeLc: Dimensional length scale (SpoungeLc_hat * Lx)
    # - sigma(x): Spatial damping coefficient (s^-1)
    # - uBC: Boundary condition state (dimensional, SI units)
    function SpoungeFun(t::Float64, x::Vector{Matrix{Float64}})
    
        # Scale spounge zone length to dimensional units
        Lc_spounge      = SpoungeLc_hat * Lx
        lambda          = abs(uR) + sqrt(gamma*RTR)
        tc_spounge      = Lc_spounge / lambda         # Characteristic wave propagation time in spounge zone
        sigma           = @tturbo @. 2.0/tc_spounge * exp(-5*(x[1]-Lx)^2/(Lc_spounge^2))
        
        # Boundary condition state (all quantities in dimensional SI units)
        uBC             = [ 0.0, 0.0, rhoR, rhoR*uR, 0.0, pR/(gamma-1.0) + 0.5*rhoR*uR^2 - rhoR*Q0 ]
        
        return sigma, uBC
        
    end
    GasModel.SpoungeBool    = SpoungeBool
    GasModel.SpoungeFun     = SpoungeFun

    # Uncomment this block to inspect the sponge profile before the simulation.
    #=
    x_sponge = reshape(collect(range(0.0, Lx, length=500)), :, 1)
    sigma_sponge, _ = SpoungeFun(0.0, [x_sponge, fill(y2/2, size(x_sponge))])
    figure()
    plot(vec(x_sponge) ./ Lx, vec(sigma_sponge) .* tc, label="sigma * tc")
    PyPlot.axvline(1.0-SpoungeLc_hat, color="gray", linestyle="--", label="Lx - sponge scale")
    xlabel("x / Lx")
    ylabel("sigma * tc")
    title(SpoungeBool ? "Sponge enabled" : "Sponge disabled (profile only)")
    legend()
    grid(true)
    PyPlot.show(block=true)
    # return  # Uncomment to stop after inspecting the profile.
    =#

    #---------------------------------------------------------------------
    #PRE-PROCESS STAGE:
    
    #Load LIRKHyp solver structure with default data. Modify the default data if necessary:
    solver                  = LIRKHyp_Start(GasModel)
    solver.ProblemName      = "ShockFlame_FXP"
    solver.SC               = SC
    solver.MeshFile         = MeshFile
    solver.nBounds          = 4
    solver.FesOrder         = FesOrder
    solver.RKMethod         = RKMethod
    solver.Deltat0          = Deltat0
    solver.tf               = tf
    solver.AMA_MaxIter      = AMA_MaxIter
    solver.AMA_SizeOrder    = AMA_SizeOrder
    solver.AMA_AnisoOrder   = AMA_AnisoOrder
    solver.TolS_max         = TolS
    solver.TolS_min         = 0.01*TolS
    solver.TolT             = TolT
#     solver.CA_max           = Inf
    solver.LSType           = "DominantILU0"
    
    #Set initial and boundary conditions:
    solver.u0fun        = FW11((x) -> u0fun(x)) 
    solver.BC           = [BCW(BC_walls), BCW(BC_outlet), BCW(BC_walls), BCW(BC_inlet)]
    
    #-----------------------------------------------------------------------------
    #INITIAL CONDITION:
    
    #Compute initial condition:
    ConvFlag            = LIRKHyp_InitialCondition!(solver, DEq_MaxIter=10)
#     CheckJacobian(solver, Plot_dQ_du=true)
#     return
    
    #Location of shock and flame and relative Mach number:
    xs, xf              = xsxf_ShockFlame2(solver; p0=p0, shock_alpha=shock_alpha)
    tv                  = [ 0.0 ]
    xsv                 = [ xs ]
    xfv                 = [ xf ]
    #Mrel and xsdot are not smoothed here:
    Mrel                = NaN
    Mrelv               = [ Mrel ]
    xsdotv              = [ NaN ]
    
    #Function to plot solution:
    figv                = Vector{Figure}(undef,3)
    if PlotFig
        figv[1]         = PyPlotSubPlots(mFig, nFig, w=wFig, h=hFig, left=0.9, right=0.4, bottom=1.1, top=1.0)
        for ii=2:length(figv)
            figv[ii]    = figure()
        end
    end    
    t_lastFig           = 0.0
    ct_SaveFig          = 0
    nb_SaveFig          = 0
    function PlotSol()
    
        ct_SaveFig      += 1
        
        if PlotFig && ( solver.t-t_lastFig>=Deltat_SaveFig || 
                        ct_SaveFig==Nt_SaveFig || solver.t==tf || solver.t==0.0 )
        
            println("Plotting")
            
            figure(figv[1].number)
    
            #Loop plot variables:
            for ii=1:length(PlotVars)
            
                PyPlot.subplot(mFig, nFig, ii)
                PyPlot.cla()

                splot_fun(x1,x2)    = @mlv x1
                v_plot              = PlotNodes(splot_fun, solver, PlotVars[ii])
                xlabel(latexstring("x_1"), fontsize=10)
                title(latexstring(LatexString(GasModel, PlotVars[ii]),
                    "; t^n=", sprintf1("%.2e", solver.t*scale_t_plot)), 
                    fontsize=10)
                    
                if MarkSF
                    plot([xf*scale_x_plot, xf*scale_x_plot], [0.9*minimum(v_plot), 1.1*maximum(v_plot)], "--r")
                    plot([xs*scale_x_plot, xs*scale_x_plot], [0.9*minimum(v_plot), 1.1*maximum(v_plot)], "--g")
                end
                
                println(PlotVars[ii], ": min=", minimum(v_plot), ", max=", maximum(v_plot))
                
            end
            if SaveFig
                savefig("$(VideosUbi)ShockFlame_FXP_SC$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
        
            figure(figv[2].number)
            PyPlot.cla()
            plot(solver.mesh)
            title(latexstring("t^n=", sprintf1("%.2e", solver.t*scale_t_plot), 
                    ", nElems=", sprintf1("%d", solver.mesh.nElems)), 
                    fontsize=10)
            if SaveFig
                savefig("$(VideosUbi)ShockFlame_FXP_Mesh$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
            
#             figure(figv[3].number)
#             PyPlot.cla()
#             semilogy(solver.tv, solver.etaSv, ".-b")
#             semilogy(solver.tv, solver.etaTv, ".-g")
#             semilogy(solver.tv, solver.etaAv, ".-r")
#             if true
#                 validv  = solver.validv .== 1
#                 semilogy(solver.tv[validv], solver.etaSv[validv], "sb")
#                 semilogy(solver.tv[validv], solver.etaTv[validv], "sg")
#                 semilogy(solver.tv[validv], solver.etaAv[validv], "sr")
#             end
#             legend(["space", "time", "algebraic"])
#             xlabel(L"t")
#             if SaveFig && solver.t==tf
#                 savefig("$(VideosUbi)ShockFlame_FXP_Errors$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
#             end
            
            figure(figv[3].number)
            PyPlot.cla()
            #
            PyPlot.subplot(1, 2, 1)
            plot(tv.*scale_t_plot, xsdotv.*scale_x_plot, "-xg")
            ylim(-0.5, 2.5)
            xlabel(L"t")
            #
            PyPlot.subplot(1, 2, 2)
            plot(tv.*scale_t_plot, Mrelv, "-xb")
            ylim(-0.5, 1.5*sqrt(gamma*RT0))
            xlabel(L"t")
            #
            if SaveFig && solver.t==tf
                savefig("$(VideosUbi)ShockFlame_FXP_Errors$(SC)_$(nb_SaveFig).png", dpi=400, pad_inches=0)
            end
            
            t_lastFig           += Deltat_SaveFig
            ct_SaveFig          = 0
            nb_SaveFig          += 1
            
            println("Plot finished")
            
        end
        return
        
    end
    PlotSol()
    
    #Function to save intermediate results:
    t_lastRes           = -Deltat_SaveRes
    ct_SaveRes          = 0
    nb_SaveRes          = 0
    function SaveSol()
        
        ct_SaveRes      += 1 
        if SaveRes && ( solver.t-t_lastRes>=Deltat_SaveRes || 
                        ct_SaveRes==Nt_SaveRes || solver.t==tf || solver.t==0.0 )
            save("$(ResUbi)LIRKHyp_SC$(SC)_$(nb_SaveRes).jld2", "StudyCase", "ShockFlame_FXP",
                "ConvFlag", ConvFlag, "solver", save(solver), 
                "epsilon", GasModel.epsilon, "nu", GasModel.nu, "kappa_rho_cv", GasModel.kappa_rho_cv, "gamma", gamma, 
                "YF0", YF0,
                "delta_flame", delta_flame, "delta_shock", delta_shock, 
                "q", q, "ML", ML, "p0", p0, "shock_alpha", shock_alpha,
                "D", GasModel.D, "BI", GasModel.BI, "BR", GasModel.BR, "BB", GasModel.BB, 
                "RTI", GasModel.RTI, "RTB", GasModel.RTB, "Q0", -GasModel.hfP, 
                "rhoCJ", rhoCJ, "pCJ", pCJ, "RTCJ", RTCJ, "uCJ", uCJ, 
                "x2", x2, "y2", y2,
                "TolS", TolS, "TolT", TolT, 
                "tv", tv, "xsv", xsv, "xfv", xfv, "Mrelv", Mrelv,
                "SpoungeLc_hat", SpoungeLc_hat)
            t_lastRes   += Deltat_SaveRes
            ct_SaveRes  = 0
            nb_SaveRes  += 1
        end
        if SaveRes && solver.t==tf
            save("$(ResUbi)LIRKHyp_SC$(SC)_info.jld2", "nb_SaveRes", nb_SaveRes-1)
        end
        return
        
    end
    SaveSol()
    
    #-----------------------------------------------------------------------------
    #MARCH IN TIME:
    
    memory_reporter = MemoryReporting.reporter()
    while solver.t<tf
    
        if uppercase(RKMethod)=="KC35"
            ConvFlag    = IRK_Step!(solver)
        else
            ConvFlag    = LIRKHyp_Step!(solver, DEq_MaxIter=0)
        end
        
        if ConvFlag<=0
            break
        end

        if memory_reporter !== nothing
            MemoryReporting.report(memory_reporter, solver.Nt)
        end
        
        #Location of shock and flame, and relative Mach number to the shock speed:
        xs, xf              = xsxf_ShockFlame2(solver; p0=p0, shock_alpha=shock_alpha)
        push!(tv, solver.t)
        push!(xsv, xs)
        push!(xfv, xf)
        Nt                  = length(tv)
        xsdot               = (xsv[Nt]-xsv[Nt-1])/(tv[Nt]-tv[Nt-1])
        Mrel                = (u0-xsdot)/sqrt(gamma*RT0)
        push!(xsdotv, xsdot)
        push!(Mrelv, Mrel)
        
        PlotSol()
        SaveSol()
        
    end
    
    return solver
    
end

#--------------------------------------------------------------------------

function pYref_ShockFlame(solver::SolverData)

    GasModel    = solver.model
    
    pv          = SolutionAtNodes(solver,GasModel,"p")
    YFv         = SolutionAtNodes(solver,GasModel,"Y_F")
    
    pref        = 0.5*(maximum(pv)+minimum(pv))
    YFref       = 0.5*(maximum(YFv)+minimum(YFv))
   
    return pref, YFref
    
end

"""
    ShockPosition_Pressure_ShockFlame2(x, p; p0=1.0, shock_alpha=0.3)

Return the first qualifying pressure-gradient peak from the left, or `nothing`.
`x` must be strictly increasing and `p` is in the same units as `p0`. This is
the production pressure criterion, independent of the diagnostic scripts.
"""
function ShockPosition_Pressure_ShockFlame2(x::AbstractVector{<:Real},
    p::AbstractVector{<:Real}; p0::Real=1.0, shock_alpha::Real=0.3)

    length(x) == length(p) || throw(ArgumentError("x and p must have the same length"))
    length(x) >= 3 || throw(ArgumentError("at least three samples are required"))
    all(isfinite, x) && all(isfinite, p) || throw(ArgumentError("x and p must be finite"))
    all(>(0), diff(x)) || throw(ArgumentError("x must be strictly increasing"))
    all(v -> isfinite(v) && v > 0, (p0, shock_alpha)) ||
        throw(ArgumentError("p0 and shock_alpha must be finite and positive"))
    threshold = p0 * (1 + shock_alpha)
    isfinite(threshold) && threshold > 0 || throw(ArgumentError("threshold is not finite and positive"))

    n = length(x)
    g = Vector{Float64}(undef, n)
    g[1] = (p[2] - p[1]) / (x[2] - x[1])
    g[end] = (p[end] - p[end-1]) / (x[end] - x[end-1])
    for i in 2:n-1
        hleft, hright = x[i] - x[i-1], x[i+1] - x[i]
        left = (p[i] - p[i-1]) / hleft
        right = (p[i+1] - p[i]) / hright
        g[i] = (hright * left + hleft * right) / (hleft + hright)
    end
    all(isfinite, g) || throw(ArgumentError("pressure gradient contains nonfinite values"))

    i = 2
    while i < n
        if g[i] > g[i-1]
            plateau_end = i
            while plateau_end < n && g[plateau_end+1] == g[i]
                plateau_end += 1
            end
            if plateau_end < n && g[plateau_end] > g[plateau_end+1]
                peak = (i + plateau_end) ÷ 2
                if g[peak] > 0 && p[peak] > threshold
                    position = Float64(x[peak])
                    if i == plateau_end
                        hl, hr = x[peak] - x[peak-1], x[peak+1] - x[peak]
                        sl = (g[peak] - g[peak-1]) / hl
                        sr = (g[peak+1] - g[peak]) / hr
                        curvature = (sr - sl) / (hl + hr)
                        if isfinite(curvature) && curvature < 0
                            vertex = x[peak] - (sl + curvature * hl) / (2 * curvature)
                            position = clamp(vertex, (x[peak-1] + x[peak])/2, (x[peak] + x[peak+1])/2)
                        end
                    end
                    return position
                end
            end
            i = plateau_end + 1
        else
            i += 1
        end
    end
    return nothing
end

function xsxf_ShockFlame2(solver::SolverData; p0::Real=1.0, shock_alpha::Real=0.3)

    GasModel    = solver.model

    #Get domain limits:
    x21         = minimum(solver.mesh.NodesCoords[:,2])
    x22         = maximum(solver.mesh.NodesCoords[:,2])
    
    #Define center line:
    x1v         = sort!(unique(solver.mesh.NodesCoords[:,1]))
    x2v         = fill(0.5*(x21+x22), length(x1v))
    
    #Interpolate solution:
    u_terp,     = SolutionCompute(solver.u, solver.fes, [x1v, x2v])
    v_terp      = DepVars(GasModel, solver.t, Vector{<:AMF64}([ x1v, x2v ]), 
                    Vector{<:AMF64}(u_terp), ["p", "Q_reac"])
    pv          = v_terp[1][1][:]
    Qv          = v_terp[2][1][:]
    
    #Get flame position (maximum heat release on the center line):
    xf          = NaN
    if any(isfinite, Qv)
        xf      = x1v[argmax(replace(Qv, NaN=>-Inf))]
    end
    position = ShockPosition_Pressure_ShockFlame2(x1v, pv; p0=p0, shock_alpha=shock_alpha)
    xs = if position === nothing
        @warn "No shock pressure-gradient peak passed the pressure threshold" t=solver.t p0 shock_alpha
        NaN
    else
        position
    end
        
    return xs, xf
    
end

function Mrel_ShockFlame(solver::SolverData, tv::Vector{Float64}, 
    xsv::Vector{Float64}, xfv::Vector{Float64})

    Nt          = length(tv)
    
    #Check if this is not the first time level:
    if Nt==1
        return NaN
    end
    
    #Check that flame has passed the shock:
    if xfv[Nt]<xsv[Nt]
        return NaN
    end
    
    #Compute shock speed:
    if Nt>2
        xsdot   = (3*xsv[Nt]-4*xsv[Nt-1]+xsv[Nt-2])/(3*tv[Nt]-4*tv[Nt-1]+tv[Nt-2])
    else
        xsdot   = (xsv[Nt]-xsv[Nt-1])/(tv[Nt]-tv[Nt-1])
    end
    
    #Compute relative Mach number at mesh nodes:
    GasModel    = solver.model
    vx          = SolutionAtNodes(solver, GasModel, "vx")
    vy          = SolutionAtNodes(solver, GasModel, "vy")
    a           = SolutionAtNodes(solver, GasModel, "a")
    Mrel        = @tturbo @. sqrt((vx-xsdot)^2 + vy^2)/a
    
    #Get Mrel at left boundary:
    x1v         = solver.mesh.NodesCoords[:,1]
    sortv       = sortperm(x1v)
    x1v         = x1v[sortv]
    Mrel        = Mrel[sortv]
    aux         = findfirst(x1v.>0.0)
    
    return Mrel[aux]
    
end
