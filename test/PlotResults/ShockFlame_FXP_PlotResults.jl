# Modular version of the plotting functions of ShockFlame_FXP_PlotResults.jl.
#
# This file is self-contained: it does NOT depend on ShockFlame_FXP_PlotResults.jl.
# It only includes the common plotting tools and the reactive-gas functions
# (FXP_Overdriven, ChapmanJouguet). The few helpers taken from the original file
# (GasFXP_ShockFlame, xdotfun, SolutionAtNodes_vxRel, SolutionAtNodes_MRel, MCJ)
# are copied here with the same names, so loading both files in the same session
# redefines them (with equivalent code).
#
# Building blocks (shared by all the plot functions):
#   LoadHistory_ShockFlame   shock/flame history of a simulation
#   LoadSnapshot_ShockFlame  everything needed to plot one saved result (SFSnapshot)
#   PlotVar_ShockFlame!      draw one variable on one axis
#   RenderSnapshot_ShockFlame! draw a list of variables on a figure (+ title)
#   ZND_ShockFlame, OverlayZND_ShockFlame!  steady ZND solution and its overlay
#
# User functions built on top of them (same names as the originals, which they replace):
#   x1Plot_ShockFlame   static / interactive (Interactive=true) plot, with
#                       optional ZND overlay (ZND=true) and per-variable log scale
#   Video_ShockFlame    movie of the saved results

include("PlotResults.jl")
include("../test_CompressibleFlow/InputData/ReactiveGasAuxFunctions.jl")


#-------------------------------------------------------------------------------
# Helpers (copied from ShockFlame_FXP_PlotResults.jl)
#-------------------------------------------------------------------------------

"Gas model of result `nb` of study case `SC`."
function GasFXP_ShockFlame(SC::Int, nb::Int)

    FileName                = GetFileName(SC, nb)

    GasModel                = GasFXP()
    GasModel.epsilon        = load(FileName, "epsilon")
    GasModel.gamma          = load(FileName, "gamma")
    GasModel.nu             = load(FileName, "nu")
    GasModel.beta           = 0.0
    GasModel.kappa_rho_cv   = load(FileName, "kappa_rho_cv")
    GasModel.hfF            = 0.0
    GasModel.hfX            = 0.0
    GasModel.hfP            = -load(FileName, "Q0")
    GasModel.D              = load(FileName, "D")
    GasModel.BR             = load(FileName, "BR")
    GasModel.BI             = load(FileName, "BI")
    GasModel.BB             = load(FileName, "BB")
    GasModel.RTI            = load(FileName, "RTI")
    GasModel.RTB            = load(FileName, "RTB")

    return GasModel

end

"""
Velocity of a front with positions `xv` at times `tv`, with the first order
formula and a minimum time step `Deltat`. The first value is NaN.
"""
function xdotfun(tv::Vector{Float64}, xv::Vector{Float64}; Deltat::Float64=1.0)

    N               = length(xv)
    xdotv           = zeros(N)
    xdotv[1]        = NaN
    for ii=2:N

        # Find the maximum jj such that tv[ii]-tv[jj] > Deltat
        jj          = ii-1
        while tv[ii]-tv[jj] <= Deltat && jj > 1
            jj      -= 1
        end

        xdotv[ii]   = (xv[ii]-xv[jj])/(tv[ii]-tv[jj])

    end

    return xdotv

end

"x1 velocity relative to the velocity `v0`, at the nodes."
function SolutionAtNodes_vxRel(solver::SolverDataSave, GasModel::GasFXP, v0::Float64)

    vx          = SolutionAtNodes(solver, GasModel, "vx")
    return @. vx-v0

end

"Mach number relative to the velocity `v0`, at the nodes."
function SolutionAtNodes_MRel(solver::SolverDataSave, GasModel::GasFXP, v0::Float64)

    vx          = SolutionAtNodes(solver, GasModel, "vx")
    vy          = SolutionAtNodes(solver, GasModel, "vy")
    a           = SolutionAtNodes(solver, GasModel, "a")
    return @. sqrt((vx-v0)^2 + vy^2)/a

end

"Chapman-Jouguet Mach number."
function MCJ(gamma::Float64, q::Float64)

    return sqrt((gamma+1)*q/2)+sqrt(1.0+(gamma+1)*q/2)

end



#-------------------------------------------------------------------------------
# Data
#-------------------------------------------------------------------------------

"Data needed to plot one saved result of a ShockFlame_FXP simulation."
struct SFSnapshot
    SC::Int
    nb::Int
    solver::Any
    GasModel::Any
    q::Float64
    ML::Float64
    t::Float64         # dimensional time
    xs::Float64         # shock position (dimensional)
    xf::Float64         # flame position (dimensional)
    xsdot::Float64      # shock velocity
    Mt::Float64         # Mach number relative to the shock
    x::Vector{Float64}  # x1 coordinate of the nodes (dimensional)
    tc::Float64         # characteristic time scale
    Lc::Float64         # characteristic length scale
    scale_x_plot::Float64  # spatial scale for plotting (physical → adimensional)
    scale_t_plot::Float64  # temporal scale for plotting (physical → adimensional)
    SpoungeLc_hat::Float64 # Spounge zone length fraction (adimensional)
    # Adimensional quantities for plotting
    t_adim::Float64        # adimensional time
    x_adim::Vector{Float64}  # adimensional x coordinates
    xs_adim::Float64       # adimensional shock position
    xf_adim::Float64       # adimensional flame position
end

"Saved result numbers for this study case, in ascending order."
function SavedResults_ShockFlame(SC::Int)
    prefix = "LIRKHyp_SC$(SC)_"
    suffix = ".jld2"
    numbers = Int[]
    for name in readdir(ResUbi)
        startswith(name, prefix) && endswith(name, suffix) || continue
        number = name[length(prefix)+1:end-length(suffix)]
        nb = tryparse(Int, number)
        nb === nothing && continue
        nb >= 0 && string(nb) == number && isfile(GetFileName(SC, nb)) && push!(numbers, nb)
    end
    sort!(numbers)
    isempty(numbers) && error("No saved ShockFlame results found for SC=$SC in $(ResUbi)")
    return numbers
end

"Return `nb`, replacing the legacy value `1000` by the latest saved result."
function ResolveNb_ShockFlame(SC::Int, nb::Int)
    nb == 1000 && return last(SavedResults_ShockFlame(SC))
    return nb
end

function FinalNb_ShockFlame(nb::Union{Int,Nothing}, alias::Union{Int,Nothing}, aliasname::String)
    if nb !== nothing && alias !== nothing && nb != alias
        throw(ArgumentError("nb and $aliasname specify different final results"))
    end
    return nb === nothing ? alias : nb
end

function ConnectKeys_ShockFlame!(fig, handler, bindings)
    owners = Dict{String,String}()
    for (name, keys) in bindings, key in keys
        isempty(key) && throw(ArgumentError("Empty key for $name"))
        haskey(owners, key) && owners[key] != name &&
            throw(ArgumentError("Key '$key' is assigned to both $(owners[key]) and $name"))
        owners[key] = name
    end
    manager = fig.canvas.manager
    if manager !== nothing && manager.key_press_handler_id !== nothing
        fig.canvas.mpl_disconnect(manager.key_press_handler_id)
        fig.canvas.mpl_connect("key_press_event", event -> begin
            if event.key !== nothing && !haskey(owners, event.key)
                PyPlot.matplotlib.backend_bases.key_press_handler(event, fig.canvas, fig.canvas.toolbar)
            end
        end)
    end
    fig.canvas.mpl_connect("key_press_event", handler)
end

CaseKeys_ShockFlame(key::String) = [lowercase(key), uppercase(key)]

"Shock and flame history (tv, xsv, xfv, xsdotv) stored in result `nb` of `SC`."
function LoadHistory_ShockFlame(SC::Int, nb::Int; Deltat::Float64=1.0)
    FileName    = GetFileName(SC, nb)
    tv          = load(FileName, "tv")
    xsv         = load(FileName, "xsv")
    xfv         = load(FileName, "xfv")
    xsdotv      = xdotfun(tv, xsv, Deltat=Deltat)
    return (tv=tv, xsv=xsv, xfv=xfv, xsdotv=xsdotv)
end

"""
    LoadSnapshot_ShockFlame(SC, nb; Deltat=1.0, history=nothing)

Load result `nb` of study case `SC`. By default the shock/flame positions and
the shock velocity are taken from the history stored in the same file. If
`history` (see `LoadHistory_ShockFlame`) is given, it is indexed with the
number of time steps stored in the file instead, which avoids recomputing the
shock velocity for every frame of a video.
"""
function LoadSnapshot_ShockFlame(SC::Int, nb::Int; Deltat::Float64=1.0, history=nothing)

    FileName    = GetFileName(SC, nb)
    GasModel    = GasFXP_ShockFlame(SC, nb)
    solver      = GetSolver(SC, nb)
    q           = load(FileName, "q")
    ML          = load(FileName, "ML")

    if history === nothing
        history = LoadHistory_ShockFlame(SC, nb, Deltat=Deltat)
        Nt      = length(history.tv)
    else
        Nt      = length(load(FileName, "tv"))
    end
    xs          = history.xsv[Nt]
    xf          = history.xfv[Nt]
    xsdot       = history.xsdotv[Nt]

    gam         = GasModel.gamma
    Mt          = (ML*sqrt(gam*(1.0+q)) - xsdot)/sqrt(gam)

    # Recalculate characteristic scales from BR (always saved)
    # Exact calculation: tc = 1/BR, Lc = uc*tc where uc=sqrt(RT0)=1.0
    tc          = 1.0 / GasModel.BR
    Lc          = tc
    scale_x_plot= 1.0 / Lc
    scale_t_plot= 1.0 / tc

    SpoungeLc_hat = jldopen(FileName, "r") do file
        haskey(file, "SpoungeLc_hat") ? file["SpoungeLc_hat"] : NaN
    end

    # Compute adimensional quantities
    x_nodes     = Vector{Float64}(solver.fes.PSpace.NodesCoords[:,1])
    t_adim      = solver.t * scale_t_plot
    x_adim      = x_nodes .* scale_x_plot
    xs_adim     = xs * scale_x_plot
    xf_adim     = xf * scale_x_plot

    return SFSnapshot(SC, nb, solver, GasModel, q, ML, solver.t, xs, xf, xsdot, Mt,
                      x_nodes, tc, Lc, scale_x_plot, scale_t_plot, SpoungeLc_hat,
                      t_adim, x_adim, xs_adim, xf_adim)

end


#-------------------------------------------------------------------------------
# Drawing
#-------------------------------------------------------------------------------

"Figure with `mFig` x `nFig` subplots and the margins used by all ShockFlame plots."
function NewFigure_ShockFlame(mFig::Int, nFig::Int; w::Float64=9.50, h::Float64=6.50)
    return PyPlotSubPlots(mFig, nFig, w=w, h=h, top=0.1, right=0.1, left=1.8,
                          bottom=0.1, tTitle=0.8, bTitle=1.0)
end

"Share the x axis of all `axes` with the first one."
function ShareX_ShockFlame!(axes)
    for i = 2:length(axes)
        axes[i].sharex(axes[1])
    end
end

"""
    Legend_ShockFlame!(ax)

Legend of the labelled curves of `ax`. Simulation curves (dots) are shown as
solid lines and ZND curves (tagged with `gid="znd"`) as dashed lines, both with
the color of the curve.
"""
function Legend_ShockFlame!(ax)

    Line2D  = PyPlot.matplotlib.lines.Line2D
    handles = Any[]
    for ln in ax.get_lines()
        lab = ln.get_label()
        (isempty(lab) || startswith(lab, "_")) && continue
        style = ln.get_gid() in ("znd", "ref") ? "--" : "-"
        push!(handles, Line2D([0], [0], color=ln.get_color(), linestyle=style,
                              linewidth=1.5, label=lab))
    end
    isempty(handles) || ax.legend(handles=handles, loc="best", fontsize=TickSize)

end

"""
    PlotVar_ShockFlame!(ax, snap, var; kwargs...)

Draw variable `var` of `snap` on axis `ax` (labels, ticks, sonic line and
shock/flame marks included). Supported values: `"Y"` (`Y_F`, `Y_X`, `Y_P`),
`"mdot"` (`-mdot_F`, `mdot_X`, `mdot_P`; note the sign of `mdot_F`), `"vxRel"`, `"uRel"` (same data as
`"vxRel"`, labelled `u'`), `"MRel"`, and any variable supported by `PlotNodes`.
Returns the plotted values (used to place the marks).

Keyword arguments:
- `MarkSF`: dashed vertical lines at the flame (red) and shock (green).
- `markersize`, `color`: marker size, and color of single-curve variables.
- `legend`: draw a legend for `"Y"` and `"mdot"` (see `Legend_ShockFlame!`).
- `labelsuffix`: LaTeX text appended to the legend labels of `"Y"` and `"mdot"`.
- `curvelabel`: legend label of single-curve variables (none if empty).
- `xlabel`: LaTeX x label.
- `useAdim`: use the characteristic length for the x coordinate and marks.
"""
function PlotVar_ShockFlame!(ax, snap::SFSnapshot, var::String;
    MarkSF::Bool=false, markersize::Real=0.2, color::String="b",
    legend::Bool=true,
    labelsuffix::String="", curvelabel::String="", xlabel::String="x_1", useAdim::Bool=false)

    solver, GasModel = snap.solver, snap.GasModel
    x = useAdim ? snap.x_adim : snap.x
    PyPlot.sca(ax)      # PlotNodes draws on the current axis
    ylab        = ""
    v_plot      = Float64[]

    if var == "vxRel" || var == "uRel"
        v_plot  = SolutionAtNodes_vxRel(solver, GasModel, snap.xsdot)
        ax.plot(x, v_plot, "."*color, markersize=markersize, label=curvelabel)
        ylab    = var == "vxRel" ? "v_x'" : "u'"
    elseif var == "MRel"
        v_plot  = SolutionAtNodes_MRel(solver, GasModel, snap.xsdot)
        ax.plot(x, v_plot, "."*color, markersize=markersize, label=curvelabel)
        ylab    = "M'"
    elseif var == "Y"
        YF      = SolutionAtNodes(solver, GasModel, "Y_F")
        YX      = SolutionAtNodes(solver, GasModel, "Y_X")
        YP      = SolutionAtNodes(solver, GasModel, "Y_P")
        ax.plot(x, YF, ".b", markersize=markersize, label=latexstring("Y_F", labelsuffix))
        ax.plot(x, YX, ".g", markersize=markersize, label=latexstring("Y_X", labelsuffix))
        ax.plot(x, YP, ".r", markersize=markersize, label=latexstring("Y_P", labelsuffix))
        ylab    = "Y"
        v_plot  = YF
    elseif var == "mdot"
        mdotF   = -SolutionAtNodes(solver, GasModel, "mdot_F")
        mdotX   = SolutionAtNodes(solver, GasModel, "mdot_X")
        mdotP   = SolutionAtNodes(solver, GasModel, "mdot_P")
        ax.plot(x, mdotF, ".r", markersize=markersize, label=latexstring("-\\dot{m}_F", labelsuffix))
        ax.plot(x, mdotX, ".g", markersize=markersize, label=latexstring("\\dot{m}_X", labelsuffix))
        ax.plot(x, mdotP, ".b", markersize=markersize, label=latexstring("\\dot{m}_P", labelsuffix))
        ylab    = "\\dot{m}"
        v_plot  = vcat(mdotF, mdotX, mdotP)
    else
        scale_x = useAdim ? snap.scale_x_plot : 1.0
        splot_fun(x1,x2) = @mlv x1 * scale_x
        v_plot  = PlotNodes(splot_fun, solver, GasModel, var, color=color, markersize=markersize)
        ylab    = LatexString(GasModel, var)
        if !isempty(curvelabel) && !isempty(ax.get_lines())
            ax.get_lines()[end].set_label(curvelabel)
        end
    end

    if legend && (var == "Y" || var == "mdot")
        Legend_ShockFlame!(ax)
    end

    ax.tick_params(axis="both", which="both", labelsize=TickSize)
    ax.set_xlabel(latexstring(useAdim ? "\\hat{x}_1" : xlabel), fontsize=10)
    ax.set_ylabel(latexstring(ylab), rotation=0, fontsize=10, labelpad=10)

    # Sonic line:
    if var == "M" || var == "MRel"
        ax.plot([minimum(x), maximum(x)], [1.0, 1.0], "--k")
    end

    # Marks are tagged so that FitY_ShockFlame! ignores them.
    if MarkSF
        isfinite(snap.xf_adim) && ax.axvline(useAdim ? snap.xf_adim : snap.xf, color="r", linestyle="--", gid="sfmark")
        isfinite(snap.xs_adim) && ax.axvline(useAdim ? snap.xs_adim : snap.xs, color="g", linestyle="--", gid="sfmark")
    end

    return v_plot

end

function RelativeAbsolute_ShockFlame(var::String)
    var == "vxRel" && return "vx"
    var == "vx" && return "vxRel"
    var == "MRel" && return "M"
    var == "M" && return "MRel"
    return var
end

"""
    FitY_ShockFlame!(ax; margin=0.05)

Set the y limits of `ax` to the range of the data inside the current x limits,
plus a relative `margin` (multiplicative on a log axis). Shock/flame marks are
ignored. Does nothing if there is no finite data in view.
"""
function FitY_ShockFlame!(ax; margin::Float64=0.05)

    xl          = ax.get_xlim()
    islog       = ax.get_yscale() == "log"
    ymin, ymax  = Inf, -Inf
    for ln in ax.get_lines()
        ln.get_gid() in ("sfmark", "ref") && continue
        xd      = Float64.(ln.get_xdata())
        yd      = Float64.(ln.get_ydata())
        for k in eachindex(xd)
            ok  = isfinite(xd[k]) && isfinite(yd[k]) && xl[1] <= xd[k] <= xl[2]
            ok && islog && (ok = yd[k] > 0.0)
            if ok
                ymin = min(ymin, yd[k])
                ymax = max(ymax, yd[k])
            end
        end
    end
    (isfinite(ymin) && isfinite(ymax)) || return

    if islog
        r   = ymax > ymin ? (ymax/ymin)^margin : 2.0
        ax.set_ylim(ymin/r, ymax*r)
    else
        pad = ymax > ymin ? margin*(ymax-ymin) : max(margin*abs(ymax), 1e-3)
        ax.set_ylim(ymin-pad, ymax+pad)
    end

end

"Title with Q, M_i, M_t and t."
function SetTitle_ShockFlame!(fig, snap::SFSnapshot; useAdim::Bool=false)
    t_display = useAdim ? snap.t_adim : snap.t
    t_label = useAdim ? "\\hat{t}" : "t"
    fig.suptitle(latexstring("\\mathcal{Q}=", snap.q, ", M_i=", snap.ML,
                             ", M_t=", sprintf1("%.2f", snap.Mt),
                             ", ", t_label, "=", sprintf1("%.2e", t_display)), fontsize=10)
end

"""
    RenderSnapshot_ShockFlame!(fig, snap, PlotVars; kwargs...)

Draw `PlotVars` of `snap` in `fig`, one per subplot of a `mFig` x `nFig` grid,
with the x axes shared, and set the title. Returns the vector of axes.

By default the figure is cleared and the axes are rebuilt. If `axes` is given,
those axes are cleared and reused (they must already share their x axis).

Keyword arguments:
- `MarkSF`, `markersize`, `PlotVarKW`: forwarded to `PlotVar_ShockFlame!`
  (`PlotVarKW` is a `NamedTuple` with extra keyword arguments).
- `PlotVarKWByVar`: per-variable keyword overrides for `PlotVarKW`.
- `logVars`: variables drawn with a logarithmic y axis.
- `xlim`: `(xmin, xmax)` applied to the shared x axis; `nothing` keeps the
  automatic limits.
- `YLims`: fixed y limits per variable, e.g. `Dict("Y"=>(1e-12, 2.0))`.
- `FitY`: fit the y axis of the variables without fixed limits to the visible data.
- `title`: set the default title.
- `overlay`: function called with the axes after all variables are drawn and
  before the limits are set (used to add the ZND curves).
"""
function RenderSnapshot_ShockFlame!(fig, snap::SFSnapshot, PlotVars::Vector{String};
    mFig::Int=2, nFig::Int=3, axes=nothing,
    MarkSF::Bool=false, markersize::Real=0.2, PlotVarKW::NamedTuple=NamedTuple(),
    PlotVarKWByVar::Dict=Dict{String,NamedTuple}(),
    logVars::Vector{String}=String[],
    xlim=nothing, YLims::Dict=Dict{String,Tuple{Float64,Float64}}(), FitY::Bool=false,
    title::Bool=true, overlay=nothing, useAdim::Bool=false)

    if axes === nothing
        PyPlot.figure(fig.number)
        fig.clf()
        axes = Any[fig.add_subplot(mFig, nFig, ii) for ii=1:length(PlotVars)]
        ShareX_ShockFlame!(axes)
    else
        for ax in axes
            ax.cla()
        end
    end

    for (ii, var) in enumerate(PlotVars)
        ax = axes[ii]
        varKW = get(PlotVarKWByVar, var, PlotVarKW)
        PlotVar_ShockFlame!(ax, snap, var, MarkSF=MarkSF, markersize=markersize, useAdim=useAdim; varKW...)
        var in logVars && ax.set_yscale("log")
    end

    overlay === nothing || overlay(axes)

    if useAdim && xlim === nothing
        xlim = (minimum(snap.x_adim), maximum(snap.x_adim))
    end
    xlim === nothing || axes[1].set_xlim(xlim[1], xlim[2])

    for (ii, var) in enumerate(PlotVars)
        if haskey(YLims, var)
            axes[ii].set_ylim(YLims[var][1], YLims[var][2])
        elseif FitY
            FitY_ShockFlame!(axes[ii])
        end
    end

    title && SetTitle_ShockFlame!(fig, snap, useAdim=useAdim)

    return axes

end


#-------------------------------------------------------------------------------
# Plot functions
#-------------------------------------------------------------------------------

"""
    x1Plot_ShockFlame(SC, nb=1000; kwargs...)

Plot `PlotVars` of result `nb` of study case `SC` (`nb=1000` is the latest
saved result, also selected when `nb` is omitted), one per subplot. The latest
result is found from the files even if the info file is absent or outdated.
This is the single user entry point replacing
`x1Plot_ShockFlame`, `x1Plot_ShockFlame_v2` and
`Compare_ShockFlame_FXP_Overdriven`.

Keyword arguments:
- `Interactive`: browse the saved results. Keys `6` / `4` go to the next /
  previous result (the arrow keys are reserved by matplotlib for the toolbar
  views) and `keyLog` (default `l`) toggles the y scale of the subplot under the
  mouse between linear and logarithmic. `keyRel` (default `r`) toggles `vxRel` /
  `vx` and `MRel` / `M` in all matching subplots. The x zoom is kept when changing
  result, the y axis is fitted to the visible data, and the Home button
  returns to the x view of the first result. Key presses received while a
  frame is being drawn are ignored. The function returns when the window is
  closed.
- `ZND`: overlay the steady ZND overdriven solution (`FXP_Overdriven`, same
  parameters, shifted to the shock position; recomputed for every result when
  combined with `Interactive`). It is only computed and drawn when the shock
  is overdriven, `M_t > M_CJ`; otherwise the simulation is plotted alone. In
  interactive mode key `keyZND` (default `z`) switches the ZND on and off
  (it can be started off with `ZND=false` and enabled later). `vxRel` is
  compared with the ZND `uRel`; absolute `vx` and `M` have no ZND counterpart.
  Variables without a ZND counterpart are plotted without it (with a warning).
- `RTcross`: draw the crossover temperature `RT_cross` (same formula as
  `compute_FXP_Tcross`, for the `M_t` of each result) as a red dashed horizontal
  line in the `RT` subplots; key `keyRTcross` (default `T`, also `t`) switches it
  on and off. It is computed for any `M_t` (with a warning if `M_t < M_CJ`).
- `MCJline`: draw `M_CJ` as a red dashed horizontal line in the Mach number
  subplots (`MRel`); key `keyMCJ` (default `M`, also `m`) switches it on and off.
- `keyRel`: interactive key for toggling relative and absolute velocity/Mach
  variables (default `r`).
- `logVars`: variables drawn with a logarithmic y axis (e.g. `["Y","p"]`).
- `FitY`: fit the y axis to the visible data. Defaults to `Interactive`.
- `useAdim`: plot space and time in characteristic units (default `true`);
  set to `false` for the stored solver coordinates and time. `Deltat` remains
  in solver time units in either mode.
- `xf_znd`, `Deltax_znd`: length and step in characteristic length units.
- `YLims`: fixed y limits per variable, e.g. `Dict("Y"=>(1e-12, 2.0))`.
- `MarkSF`, `markersize`, `mFig`, `nFig`, `w`, `h`, `Deltat`: as in the
  original functions.
- `SaveFig`: save `\$(FigUbi)SC<SC>_<nb>_x1Plot.png` (`..._x1Plot_ZND.png` with
  `ZND=true`) for the initial result and when navigating to another result.
  Interactive display toggles do not trigger another high-resolution save.

See `PlotVar_ShockFlame!` and `RenderSnapshot_ShockFlame!` for details.
Returns `nothing`.
"""
function x1Plot_ShockFlame(SC::Int, nb::Int=1000;
    PlotVars::Vector{String}=["Y", "p", "vxRel", "mdot", "RT", "MRel"],
    MarkSF::Bool=false,
    mFig::Int=2, nFig::Int=3,
    SaveFig::Bool=false, w::Float64=9.50, h::Float64=6.50,
    Deltat::Float64=1.0,
    logVars::Vector{String}=String[], FitY::Union{Bool,Nothing}=nothing,
    YLims::Dict=Dict{String,Tuple{Float64,Float64}}(), markersize::Real=0.2,
    Interactive::Bool=false,
    ZND::Bool=false, xf_znd::Float64=20.0, Deltax_znd::Float64=1e-3,
    keysNext::Vector{String}=["6"], keysPrev::Vector{String}=["4"], keyLog::String="l", keyZND::String="z",
    RTcross::Bool=false, keyRTcross::String="T", MCJline::Bool=false, keyMCJ::String="M", keyRel::String="r",
    useAdim::Bool=true)

    ZND && CheckZNDVars_ShockFlame(PlotVars)
    activePlotVars = copy(PlotVars)
    zndOn        = Ref(ZND)
    rtcOn        = Ref(RTcross)
    mcjOn        = Ref(MCJline)
    fit          = FitY === nothing ? Interactive : FitY
    logSet       = Set(logVars)
    nb_ref       = Ref{Int}(ResolveNb_ShockFlame(SC, nb))
    first_render = Ref(true)
    busy         = Ref(false)
    home_xlims   = Ref(Any[])
    current_snap = Ref{Union{Nothing,SFSnapshot}}(nothing)
    plot_axes    = Ref{Any}(nothing)
    znd_cache    = Ref{Any}(nothing)
    znd_cached   = Ref(false)
    rtc_cache    = Ref{Any}(nothing)
    rtc_cached   = Ref(false)
    fig          = NewFigure_ShockFlame(mFig, nFig, w=w, h=h)

    function render_plot(current_nb::Int; save_frame::Bool=true)

        if !isfile(GetFileName(SC, current_nb))
            @warn "File does not exist: $(GetFileName(SC, current_nb))"
            return false
        end

        # x view to keep (only if a previous render exists).
        prev_xlims  = first_render[] ? nothing : [ax.get_xlim() for ax in plot_axes[]]

        snap = current_snap[]
        if snap === nothing || snap.nb != current_nb
            snap = LoadSnapshot_ShockFlame(SC, current_nb, Deltat=Deltat)
            current_snap[] = snap
            znd_cache[] = nothing
            znd_cached[] = false
            rtc_cache[] = nothing
            rtc_cached[] = false
        end

        # The ZND solution is only computed and drawn for overdriven shocks (M_t > M_CJ).
        mcj         = MCJ(snap.GasModel.gamma, snap.q)
        hasZNDVars  = any(var -> ZNDName_ShockFlame(var) !== nothing, activePlotVars)
        useZND      = zndOn[] && hasZNDVars && snap.Mt > mcj
        if zndOn[] && hasZNDVars && snap.Mt <= mcj
            @info "ZND not drawn for nb=$(current_nb): M_t=$(round(snap.Mt, digits=3)), M_CJ=$(round(mcj, digits=3))"
        end
        if useZND && !znd_cached[]
            znd_cache[] = ZND_ShockFlame(snap, xf_znd=xf_znd, Deltax_znd=Deltax_znd)
            znd_cached[] = true
        end
        znd         = useZND ? znd_cache[] : nothing
        dynKW       = (color="r", legend=false, labelsuffix="~(\\mathrm{dyn})",
                       curvelabel="dyn", xlabel="x")
        dynKWByVar  = useZND ? Dict(var => dynKW for var in activePlotVars
                                    if ZNDName_ShockFlame(var) !== nothing) : Dict{String,NamedTuple}()
        if rtcOn[] && !rtc_cached[]
            rtc_cache[] = RTCross_ShockFlame(snap)
            rtc_cached[] = true
        end
        rtc         = rtcOn[] ? rtc_cache[] : nothing
        overlay     = (useZND || rtc !== nothing || mcjOn[]) ? (axs -> begin
                        useZND && OverlayZND_ShockFlame!(axs, activePlotVars, znd, scale_x=useAdim ? snap.scale_x_plot : 1.0)
                        rtc === nothing || OverlayRTCross_ShockFlame!(axs, activePlotVars, rtc)
                        mcjOn[] && OverlayMCJ_ShockFlame!(axs, activePlotVars, mcj)
                      end) : nothing
        axes        = RenderSnapshot_ShockFlame!(fig, snap, activePlotVars, mFig=mFig, nFig=nFig,
                        axes=plot_axes[], PlotVarKWByVar=dynKWByVar,
                        MarkSF=MarkSF, markersize=markersize, logVars=collect(logSet),
                        YLims=YLims, FitY=false,
                        overlay=overlay, useAdim=useAdim)
        plot_axes[] = axes

        # The toolbar view stack is rebuilt: [home view, current view].
        first_render[] && (home_xlims[] = [ax.get_xlim() for ax in axes])
        toolbar     = fig.canvas.toolbar
        if toolbar !== nothing
            toolbar.update()
            toolbar.push_current()
        end
        if prev_xlims !== nothing
            for (idx, ax) in enumerate(axes)
                ax.set_xlim(prev_xlims[idx])
            end
        end
        if fit
            for (idx, ax) in enumerate(axes)
                haskey(YLims, activePlotVars[idx]) || FitY_ShockFlame!(ax)
            end
        end
        toolbar !== nothing && toolbar.push_current()

        if SaveFig && save_frame
            fig.savefig("$(FigUbi)SC$(SC)_$(current_nb)_x1Plot$(useZND ? "_ZND" : "").png",
                        dpi=800, pad_inches=0)
        end

        first_render[] = false
        fig.canvas.draw_idle()
        return true

    end

    function on_key(event)
        # Ignore key presses while a frame is being drawn, so they do not pile up.
        busy[] && return
        key     = event.key
        delta   = key in keysNext ? 1 : (key in keysPrev ? -1 : 0)
        target  = nb_ref[] + delta
        toggle  = key == keyLog
        toggleZ = key == keyZND
        toggleT = key !== nothing && lowercase(key) == lowercase(keyRTcross)
        toggleM = key !== nothing && lowercase(key) == lowercase(keyMCJ)
        toggleR = key == keyRel
        if toggle || toggleZ || toggleT || toggleM || toggleR || (delta != 0 && target >= 1)
            busy[] = true
            try
                if toggleM
                    mcjOn[] = !mcjOn[]
                    render_plot(nb_ref[]; save_frame=false)
                elseif toggleT
                    rtcOn[] = !rtcOn[]
                    render_plot(nb_ref[]; save_frame=false)
                elseif toggleZ
                    zndOn[] = !zndOn[]
                    zndOn[] && CheckZNDVars_ShockFlame(PlotVars)
                    render_plot(nb_ref[]; save_frame=false)
                elseif toggle
                    # Toggle the y scale without rebuilding the plotted data.
                    inax    = event.inaxes
                    idx     = inax === nothing ? nothing : findfirst(a -> a == inax, plot_axes[])
                    if idx !== nothing && idx <= length(activePlotVars)
                        var = activePlotVars[idx]
                        var in logSet ? delete!(logSet, var) : push!(logSet, var)
                        ax = plot_axes[][idx]
                        ax.set_yscale(var in logSet ? "log" : "linear")
                        if haskey(YLims, var)
                            ax.set_ylim(YLims[var][1], YLims[var][2])
                        elseif fit
                            FitY_ShockFlame!(ax)
                        end
                        fig.canvas.draw_idle()
                    end
                elseif toggleR
                    activePlotVars .= RelativeAbsolute_ShockFlame.(activePlotVars)
                    render_plot(nb_ref[]; save_frame=false)
                elseif render_plot(target)
                    nb_ref[] = target
                    println("Step updated: nb = ", nb_ref[])
                end
            finally
                busy[] = false
            end
        end
    end

    render_plot(nb_ref[])
    Interactive || return nothing

    ConnectKeys_ShockFlame!(fig, on_key, [
        "next" => keysNext, "previous" => keysPrev, "log" => [keyLog],
        "ZND" => [keyZND], "RTcross" => CaseKeys_ShockFlame(keyRTcross),
        "MCJ" => CaseKeys_ShockFlame(keyMCJ), "relative" => [keyRel]])

    # Keep processing GUI events until the window is closed.
    while PyPlot.plt.fignum_exists(fig.number)
        pause(0.05)
    end

    return nothing

end

"""
    ZND_ShockFlame(snap; xf_znd=20.0, Deltax_znd=1e-4)

Steady ZND overdriven solution (`FXP_Overdriven`) with the same thermodynamic
and kinetic parameters as `snap`, and transmitted shock Mach number computed
from the shock velocity of the simulation. The solution is shifted to the shock
position. Returns a `NamedTuple` with the abscissa `x` and the `curves`, a
`Dict` with the ZND values of `p`, `rho`, `RT`, `uRel`, `MRel`, `YF`, `YX`, `YP`,
`mdotF` (that is, `-mdot_F`), `mdotX` and `mdotP`.
"""
function ZND_ShockFlame(snap::SFSnapshot; xf_znd::Float64=20.0, Deltax_znd::Float64=1e-4)

    FileName    = GetFileName(snap.SC, snap.nb)
    GasModel    = snap.GasModel
    q           = snap.q
    YF0         = load(FileName, "YF0")
    gamma       = GasModel.gamma

    # Same CJ reference state as ShockFlame_FXP and FXP_Overdriven.
    _, _, _, _, rhoCJ, _, pCJ, _ = ChapmanJouguet(gamma, q)
    RTCJ        = pCJ/rhoCJ
    RThatI      = GasModel.RTI/RTCJ
    RThatB      = GasModel.RTB/RTCJ
    BIhat       = GasModel.BI/GasModel.BR
    BBhat       = GasModel.BB*rhoCJ/GasModel.BR

    # Transmitted shock Mach number, from the shock velocity of the simulation.
    MRel        = SolutionAtNodes_MRel(snap.solver, GasModel, snap.xsdot)
    xv_znd, rhov_znd, uv_znd, pv_znd, YFv_znd, YXv_znd = FXP_Overdriven_Adaptive(
        gamma=gamma, q=q, M0=MRel[1],
        RThatI=RThatI, RThatB=RThatB, BIhat=BIhat, BBhat=BBhat,
        Deltax=Deltax_znd, xf=xf_znd, YF0=YF0, YX0=1.0-YF0)

    xv_znd      = xv_znd .* snap.Lc .+ snap.xs
    RTv_znd     = pv_znd ./ rhov_znd
    Mv_znd      = uv_znd ./ sqrt.(gamma .* RTv_znd)
    YPv_znd     = 1.0 .- YFv_znd .- YXv_znd

    kI_znd      = @. GasModel.BI*exp(-GasModel.RTI/RTv_znd)
    kB_znd      = @. GasModel.BB*exp(-GasModel.RTB/RTv_znd)
    rhoYF_znd   = rhov_znd .* YFv_znd
    rhoYX_znd   = rhov_znd .* YXv_znd
    mdotFv_znd  = @. rhoYF_znd*kI_znd + rhoYF_znd*rhoYX_znd*kB_znd     # -mdot_F
    mdotPv_znd  = @. GasModel.BR*rhoYX_znd
    mdotXv_znd  = @. mdotFv_znd - mdotPv_znd

    curves      = Dict("p"=>pv_znd, "rho"=>rhov_znd, "RT"=>RTv_znd,
                       "uRel"=>uv_znd, "MRel"=>Mv_znd,
                       "YF"=>YFv_znd, "YX"=>YXv_znd, "YP"=>YPv_znd,
                       "mdotF"=>mdotFv_znd, "mdotX"=>mdotXv_znd, "mdotP"=>mdotPv_znd)

    return (x=xv_znd, curves=curves)

end

"""
Crossover temperature `RT_cross` for the shock Mach number `M_t` of `snap`, using
the gas model parameters of the simulation. Same formula as `compute_FXP_Tcross`
but evaluated for any `M_t` (a warning is shown if `M_t < M_CJ`). Returns
`nothing` (with a warning) only if the branching rate is too low for the
formula to be defined.
"""
function RTCross_ShockFlame(snap::SFSnapshot)

    GasModel    = snap.GasModel
    gamma       = GasModel.gamma
    YF0         = load(GetFileName(snap.SC, snap.nb), "YF0")
    M0CJ, rho0, _, _, rhoCJ, _, pCJ, _ = ChapmanJouguet(gamma, snap.q)
    RTCJ        = pCJ/rhoCJ
    RThatB      = GasModel.RTB/RTCJ
    BBhat       = GasModel.BB*rhoCJ/GasModel.BR
    M0          = Float64(SolutionAtNodes_MRel(snap.solver, GasModel, snap.xsdot)[1])

    M0 < M0CJ && @warn "RT_cross computed for nb=$(snap.nb) with M_t=$(round(M0, digits=3)) < M_CJ=$(round(M0CJ, digits=3))"

    # Post-shock density for the transmitted Mach number M0 (Rankine-Hugoniot).
    rho1        = rho0*((gamma+1.0)*M0^2)/(2.0+(gamma-1.0)*M0^2)
    arg         = BBhat*YF0*rho1/rhoCJ
    if !(arg > 1.0)
        @warn "RT_cross not drawn for nb=$(snap.nb): BBhat too low to sustain chain-branching"
        return nothing
    end
    return RThatB*RTCJ/log(arg)

end

"Draw `M_CJ` as a red dashed horizontal line in the Mach number (`MRel`) subplots."
function OverlayMCJ_ShockFlame!(axes, PlotVars::Vector{String}, MCJval::Real)

    for (ii, var) in enumerate(PlotVars)
        var in ("MRel", "M") || continue
        axes[ii].axhline(MCJval, color="r", linestyle="--", gid="ref",
                         label=latexstring("M_{CJ}"))
        Legend_ShockFlame!(axes[ii])
    end

end

"Draw `RT_cross` as a red dashed horizontal line in the `RT` subplots."
function OverlayRTCross_ShockFlame!(axes, PlotVars::Vector{String}, RTc::Real)

    for (ii, var) in enumerate(PlotVars)
        var == "RT" || continue
        axes[ii].axhline(RTc, color="r", linestyle="--", gid="ref",
                         label=latexstring("RT_{cross}"))
        Legend_ShockFlame!(axes[ii])
    end

end

"ZND curve name of variable `var`, or `nothing` if it has none (`vxRel` is compared with `uRel`)."
function ZNDName_ShockFlame(var::String)
    var == "vxRel" && return "uRel"
    return var in ["Y", "mdot", "p", "rho", "RT", "uRel", "MRel"] ? var : nothing
end

"Warn about the variables of `PlotVars` that have no ZND counterpart (drawn without it)."
function CheckZNDVars_ShockFlame(PlotVars::Vector{String})
    for vble in PlotVars
        ZNDName_ShockFlame(vble) === nothing &&
            @warn "Variable $(vble) has no ZND counterpart: plotted without it"
    end
end

"Add the ZND curves `znd` (see `ZND_ShockFlame`) to the axes drawn for `PlotVars`."
function OverlayZND_ShockFlame!(axes, PlotVars::Vector{String}, znd; scale_x::Float64=1.0)

    c   = znd.curves
    x   = znd.x .* scale_x
    for (ii, var) in enumerate(PlotVars)
        ax = axes[ii]
        if var == "Y"
            ax.plot(x, c["YF"], "--b", gid="znd", label=latexstring("Y_F~(\\mathrm{znd})"))
            ax.plot(x, c["YX"], "--g", gid="znd", label=latexstring("Y_X~(\\mathrm{znd})"))
            ax.plot(x, c["YP"], "--r", gid="znd", label=latexstring("Y_P~(\\mathrm{znd})"))
        elseif var == "mdot"
            ax.plot(x, c["mdotF"], "--r", gid="znd", label=latexstring("-\\dot{m}_F~(\\mathrm{znd})"))
            ax.plot(x, c["mdotX"], "--g", gid="znd", label=latexstring("\\dot{m}_X~(\\mathrm{znd})"))
            ax.plot(x, c["mdotP"], "--b", gid="znd", label=latexstring("\\dot{m}_P~(\\mathrm{znd})"))
        else
            name = ZNDName_ShockFlame(var)
            name === nothing && continue
            # Same color as the simulation curve of this subplot (red, see dynKW).
            ax.plot(x, c[name], "--r", gid="znd", label="znd")
        end
        Legend_ShockFlame!(ax)
    end

end

"""
    Video_ShockFlame(SC; kwargs...)

Plays saved results from `n1` through `nb` (inclusive), or through the latest
saved result when `nb` is omitted. `n2` is a backward-compatible alias for `nb`;
if both are provided, they must agree. Missing intermediate files are skipped
with a warning; an explicitly requested final result must exist. Discovery does not require
the info file. Space or `p` pauses/resumes, `keyRestart` (default `i`) restarts,
and `keyRel` (default `r`) toggles relative/absolute velocity and Mach subplots
as in `x1Plot`. When front data are unavailable a frame is still displayed, but
relative variables or tracking limits may be unavailable. When the video ends
the last frame stays until restarted or the window is closed.

Extra keyword arguments with respect to the original: `logVars` (default
empty: all subplots start in linear scale), `YLims` (default `Y` in `[1e-12, 2]`,
applied only while that subplot is in log scale) and `FitY` (fit the y axis
of the other subplots to the visible data, default `true`). `TrackSF` keeps
the x axis between the shock and the flame. Key `keyLog` (default `l`) toggles
the y scale of the subplot under the mouse between linear and logarithmic; the
choice is kept in the following frames. `RTcross` (default `false`) draws `RT_cross`
in the `RT` subplots; key `keyRTcross` (default `T`) switches it on and off.
`MCJline` (default `false`) draws `M_CJ` in the `MRel` subplots; key `keyMCJ`
(default `M`) switches it on and off. `useAdim` defaults to `true` and displays
space and time in characteristic units; `false` uses the saved solver units.
`Deltat` is always in solver time units. The `TrackSF` margin is 0.2 in solver
length units, converted along with the front positions.
"""
function Video_ShockFlame(SC::Int;
    PlotVars::Vector{String}=["Y", "p", "vxRel", "mdot", "RT", "MRel"],
    MarkSF::Bool=false, TrackSF::Bool=true,
    mFig::Int=2, nFig::Int=3,
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50,
    tPlot::Float64=30.0, n1::Int=0, nb::Union{Int,Nothing}=nothing,
    n2::Union{Int,Nothing}=nothing, Deltat::Float64=1.0,
    logVars::Vector{String}=String[], keyLog::String="l", keyRel::String="r", keyRestart::String="i",
    RTcross::Bool=false, keyRTcross::String="T", MCJline::Bool=false, keyMCJ::String="M",
    YLims::Dict=Dict("Y"=>(1e-12, 2.0)), FitY::Bool=true, markersize::Real=0.2,
    useAdim::Bool=true)

    n1 >= 0 || throw(ArgumentError("n1 must be nonnegative"))
    final_nb    = FinalNb_ShockFlame(nb, n2, "n2")
    saved       = SavedResults_ShockFlame(SC)
    nbMax       = final_nb === nothing ? last(saved) : final_nb
    nbMax >= n1 || throw(ArgumentError("nb must be at least n1"))
    nbMax in saved || throw(ArgumentError("Saved result nb=$nbMax does not exist for SC=$SC"))
    frames      = filter(nb -> n1 <= nb <= nbMax, saved)
    missing     = setdiff(n1:nbMax, frames)
    isempty(missing) || @warn "Video_ShockFlame: missing saved results for SC=$SC; skipping" missing
    history     = LoadHistory_ShockFlame(SC, nbMax, Deltat=Deltat)
    tPlot_nFig  = tPlot/length(frames)

    fig         = NewFigure_ShockFlame(mFig, nFig, w=w, h=h)
    axes        = fig.axes[1:length(PlotVars)]
    ShareX_ShockFlame!(axes)
    activePlotVars = copy(PlotVars)

    is_paused    = Ref(false)
    restart_flag = Ref(false)
    logSet       = Set(logVars)
    last_snap    = Ref{Any}(nothing)
    rtcOn        = Ref(RTcross)
    mcjOn        = Ref(MCJline)

    function draw_frame(snap)
        rtc     = rtcOn[] ? RTCross_ShockFlame(snap) : nothing
        scale_x = useAdim ? snap.scale_x_plot : 1.0
        track_xlim = TrackSF && isfinite(snap.xs) && isfinite(snap.xf)
        RenderSnapshot_ShockFlame!(fig, snap, activePlotVars, axes=axes,
            MarkSF=MarkSF, markersize=markersize, logVars=collect(logSet),
            xlim=track_xlim ? ((snap.xs-0.2)*scale_x, (snap.xf+0.2)*scale_x) : nothing,
            YLims=filter(kv -> kv.first in logSet, YLims), FitY=FitY,
            overlay=(rtc === nothing && !mcjOn[]) ? nothing : (axs -> begin
                rtc === nothing || OverlayRTCross_ShockFlame!(axs, activePlotVars, rtc)
                mcjOn[] && OverlayMCJ_ShockFlame!(axs, activePlotVars, MCJ(snap.GasModel.gamma, snap.q))
            end), useAdim=useAdim)
    end

    function on_key(event)
        if event.key == " " || event.key == "p"
            is_paused[] = !is_paused[]
        elseif event.key !== nothing && lowercase(event.key) == lowercase(keyMCJ)
            mcjOn[] = !mcjOn[]
            if last_snap[] !== nothing
                draw_frame(last_snap[])
                fig.canvas.draw_idle()
            end
        elseif event.key !== nothing && lowercase(event.key) == lowercase(keyRTcross)
            rtcOn[] = !rtcOn[]
            if last_snap[] !== nothing
                draw_frame(last_snap[])
                fig.canvas.draw_idle()
            end
        elseif event.key == keyRestart
            restart_flag[] = true
            is_paused[] = false
        elseif event.key == keyRel
            activePlotVars .= RelativeAbsolute_ShockFlame.(activePlotVars)
            if last_snap[] !== nothing
                draw_frame(last_snap[])
                fig.canvas.draw_idle()
            end
        elseif event.key == keyLog
            # Toggle the y scale of the subplot under the mouse; it is kept in
            # the following frames, and the current frame is redrawn.
            inax    = event.inaxes
            idx     = inax === nothing ? nothing : findfirst(a -> a == inax, axes)
            if idx !== nothing
                var = activePlotVars[idx]
                var in logSet ? delete!(logSet, var) : push!(logSet, var)
                if last_snap[] !== nothing
                    draw_frame(last_snap[])
                    fig.canvas.draw_idle()
                end
            end
        end
    end
    ConnectKeys_ShockFlame!(fig, on_key, [
        "pause" => [" ", "p"], "restart" => [keyRestart], "relative" => [keyRel],
        "log" => [keyLog], "RTcross" => CaseKeys_ShockFlame(keyRTcross),
        "MCJ" => CaseKeys_ShockFlame(keyMCJ)])

    kSave       = 0
    while true
        restart_flag[] = false

        for nb in frames

            PyPlot.plt.fignum_exists(fig.number) || return
            t0      = time()
            snap    = LoadSnapshot_ShockFlame(SC, nb, history=history)

            if !(isfinite(snap.xs) && isfinite(snap.xf) && isfinite(snap.xsdot))
                @warn "Video_ShockFlame: nb=$nb has unavailable front data; relative values or tracking limits may be unavailable" xs=snap.xs xf=snap.xf xsdot=snap.xsdot
            end

            last_snap[] = snap
            draw_frame(snap)

            if SaveFig
                fig.savefig("$(FigUbi)SC$(SC)_$(kSave).png", dpi=800, pad_inches=0)
            end
            kSave   += 1

            pause(min(max(0.01, tPlot_nFig-(time()-t0)), 0.5))

            while is_paused[]
                pause(0.05)
                PyPlot.plt.fignum_exists(fig.number) || return
            end

            restart_flag[] && break
        end

        # Keep the last frame until restarted or the window is closed.
        while !restart_flag[]
            PyPlot.plt.fignum_exists(fig.number) || break
            pause(0.1)
        end

        restart_flag[] || break
    end

    return nothing

end


function Contour_ShockFlame(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50, 
    PlotVars::Vector{String}=["Y_F", "v", "p", "RT", "rho", "M"], mFig::Int=3, nFig::Int=2, 
    PlotTitle::Bool=false)

    GasModel    = GasFXP_ShockFlame(SC, nb)
    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    
    for ii=1:length(PlotVars)
        
        if PlotTitle
            PyPlotFigure(w=w, h=h+0.8, top=0.8, bottom=0.0, left=0.0, right=0.0)
        else
            PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
        end
    
        #Numerical solution:
        PlotContour(solver, GasModel, PlotVars[ii], delta=1e-5)
        tick_params(axis="both", which="both", labelsize=TickSize)
        axis("off")
        if PlotTitle
            title(latexstring(LatexString(GasModel, PlotVars[ii]),
                "; t^n=", sprintf1("%.2e", solver.t)), 
                fontsize=10)
        end
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_$(nb)_Contour_$(PlotVars[ii]).png", dpi=800, pad_inches=0)
        end
    
    end
    
    return
    
end

function Mesh_ShockFlame(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)


    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    
    PyPlotFigure(w=w, h=h, bottom=1.5)
    PlotContour(solver.u[1], solver.fes)
    PlotMesh!(SC, nb, color="w")
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Mesh.png", dpi=800, pad_inches=0)
    end
    
    return
    
end


"""
    Track_ShockFlame(SC; n1=nothing, nb=nothing, useAdim=true, Deltat=1.0, ...)

Plot front positions, velocities, separation and Mach number over time.
The history comes from `nb`, or from the latest saved result when omitted.
`n2` and `nbMax` are aliases for `nb` (`nbMax < 0` is ignored); supplied final
indices must agree. If `n1` is given, plot from that saved result's time,
inclusive. If its file is missing, use the next available result in the range
with a warning. Without `n1`, plot the complete history. Velocities are computed
from the complete history before restricting the displayed interval. The
separate info file is not required.
By default space and time use the characteristic scales reconstructed from
`BR`; `useAdim=false` displays solver units. Velocities are divided by the
characteristic velocity. `Deltat` is always specified in solver time units.
"""
function Track_ShockFlame(SC::Int; SaveFig::Bool=false, w::Float64=10.0, h::Float64=10.0,
    Deltat::Float64=1.0, n1::Union{Int,Nothing}=nothing, nb::Union{Int,Nothing}=nothing,
    n2::Union{Int,Nothing}=nothing, nbMax::Int=-1, useAdim::Bool=true)

    nb !== nothing && nb < 0 && throw(ArgumentError("nb must be nonnegative"))
    n2 !== nothing && n2 < 0 && throw(ArgumentError("n2 must be nonnegative"))
    n1 !== nothing && n1 < 0 && throw(ArgumentError("n1 must be nonnegative"))
    final_nb = FinalNb_ShockFlame(nb, n2, "n2")
    final_nb = FinalNb_ShockFlame(final_nb, nbMax < 0 ? nothing : nbMax, "nbMax")
    nbMax = final_nb === nothing ? last(SavedResults_ShockFlame(SC)) : final_nb
    n1 !== nothing && n1 > nbMax && throw(ArgumentError("n1 must not exceed the final result nb=$nbMax"))
    
    FileName    = GetFileName(SC, nbMax)
    q           = load(FileName, "q")
    ML          = load(FileName, "ML")
    tv          = load(FileName, "tv")
    xsv         = load(FileName, "xsv")
    xfv         = load(FileName, "xfv")
    gamma       = load(FileName, "gamma")
    BR          = load(FileName, "BR")
    tc          = 1.0 / BR
    Lc          = tc
    scale_x     = useAdim ? 1.0 / Lc : 1.0
    scale_t     = useAdim ? 1.0 / tc : 1.0
    t_label     = useAdim ? "\\hat{t}" : "t"
#     Mtv         = load(FileName, "Mrelv")

    #Relative Mach number:
    u0              = ML*sqrt(gamma*(1.0+q))    #u0=uL, RT0=1, RTL=RT0*(1+q)
    a0              = sqrt(gamma*1.0)
    xsdotv          = xdotfun(tv, xsv, Deltat=Deltat)
    xfdotv          = xdotfun(tv, xfv, Deltat=Deltat)
    Mtv             = @. (u0-xsdotv)/a0
#     Mtv[xfv.<xsv]   .= NaN

    if n1 !== nothing
        saved = SavedResults_ShockFlame(SC)
        first_idx = findfirst(index -> n1 <= index <= nbMax, saved)
        first_idx === nothing && throw(ArgumentError("No saved results between n1=$n1 and nb=$nbMax"))
        start_nb = saved[first_idx]
        start_nb == n1 || @warn "Track_ShockFlame: result n1=$n1 is missing; starting at nb=$start_nb"
        start_tv = load(GetFileName(SC, start_nb), "tv")
        start_index = length(start_tv)
        1 <= start_index <= length(tv) && isequal(start_tv[end], tv[start_index]) ||
            error("Track_ShockFlame: inconsistent time histories for nb=$start_nb and nb=$nbMax")
        range = start_index:length(tv)
        tv, xsv, xfv = tv[range], xsv[range], xfv[range]
        xsdotv, xfdotv, Mtv = xsdotv[range], xfdotv[range], Mtv[range]
    end
    tv_plot = tv .* scale_t
    
    #Plot xs(t) and xf(t):
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.0, left=1.5, right=0.2)
    plot(tv_plot, xsv .* scale_x, "b")
    plot(tv_plot, xfv .* scale_x, "r")
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring(t_label), fontsize=10)
    ylabel(latexstring(useAdim ? "\\hat{x}_1" : "x_1"), fontsize=10, rotation=0)
    legend(["shock", "flame"], fontsize=8)
#     title(latexstring("\\mathcal{Q}=", q, ", \\mathcal{M}_i=", ML), fontsize=10)
    title(latexstring("\\mathcal{Q}=", q, ", M_L=", ML), fontsize=10)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_TrackShockFlame.png", dpi=800, pad_inches=0)
    end
    
    #Plot xsdot:
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.0, left=1.5, right=0.2)
    xsdot_plot      = xsdotv .* (scale_x / scale_t)
    xfdot_plot      = xfdotv .* (scale_x / scale_t)
    plot(tv_plot, xsdot_plot, "b")
    plot(tv_plot, xfdot_plot, "r")
    legend([latexstring(useAdim ? "\\dot{\\hat{x}}_s" : "\\dot{x}_s"),
            latexstring(useAdim ? "\\dot{\\hat{x}}_f" : "\\dot{x}_f")], fontsize=8)
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring(t_label), fontsize=10)
    ylabel(latexstring(useAdim ? "\\dot{\\hat{x}}" : "\\dot{x}"), fontsize=10, rotation=0)
    allv            = vcat(xsdot_plot, xfdot_plot)
    finite_velocities = allv[isfinite.(allv)]
    isempty(finite_velocities) || ylim(-0.1, 1.2*maximum(finite_velocities))
    title(latexstring("\\mathcal{Q}=", q, ", \\mathcal{M}_i=", ML), fontsize=10)
    grid("on")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_TrackShockFlame2.png", dpi=800, pad_inches=0)
    end
    
    #Plot xs(t) and xf(t):
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.0, left=1.5, right=0.2)
    auxv            = xfv-xsv
    auxv[xfv.<xsv]  .= NaN
    plot(tv_plot, auxv .* scale_x, "m")
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring(t_label), fontsize=10)
    ylabel(latexstring(useAdim ? "\\hat{x}_f-\\hat{x}_s" : "x_f-x_s"), fontsize=10, rotation=0)
    title(latexstring("\\mathcal{Q}=", q, ", \\mathcal{M}_i=", ML), fontsize=10)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_TrackShockFlame3.png", dpi=800, pad_inches=0)
    end
    
    #Plot Mrel(t):
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.0, left=1.5, right=0.2)
    plot(tv_plot, Mtv, "b")
    plot([tv_plot[1], tv_plot[end]], fill(MCJ(gamma, q), 2), "--r")
    legend([latexstring("\\mathcal{M}_t"), latexstring("\\mathcal{M}_{CJ}")], fontsize=8)
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring(t_label), fontsize=10)
    ylabel(latexstring("\\mathcal{M}"), fontsize=10, rotation=0, labelpad=10.0)    
    title(latexstring("\\mathcal{Q}=", q, ", \\mathcal{M}_i=", ML), fontsize=10)
    finite_mach = Mtv[isfinite.(Mtv)]
    ylim(-0.1, 1.2*(isempty(finite_mach) ? MCJ(gamma, q) : max(MCJ(gamma, q), maximum(finite_mach))))
    grid("on")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Mt.png", dpi=800, pad_inches=0)
    end
#     println("M_rel=", Mtv)
    
    return
    
end




