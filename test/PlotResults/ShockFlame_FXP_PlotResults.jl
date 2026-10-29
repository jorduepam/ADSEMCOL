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
    t::Float64
    xs::Float64         # shock position
    xf::Float64         # flame position
    xsdot::Float64      # shock velocity
    Mt::Float64         # Mach number relative to the shock
    x::Vector{Float64}  # x1 coordinate of the nodes
end

"Return `nb`, replacing `nb==1000` by the latest saved result of study case `SC`."
function ResolveNb_ShockFlame(SC::Int, nb::Int)
    nb == 1000 && return load("$(ResUbi)LIRKHyp_SC$(SC)_info.jld2", "nb_SaveRes")
    return nb
end

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

    return SFSnapshot(SC, nb, solver, GasModel, q, ML, solver.t, xs, xf, xsdot, Mt,
                      Vector{Float64}(solver.fes.PSpace.NodesCoords[:,1]))

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
"""
function PlotVar_ShockFlame!(ax, snap::SFSnapshot, var::String;
    MarkSF::Bool=false, markersize::Real=0.2, color::String="b",
    legend::Bool=true,
    labelsuffix::String="", curvelabel::String="", xlabel::String="x_1")

    solver, GasModel, x = snap.solver, snap.GasModel, snap.x
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
        splot_fun(x1,x2) = @mlv x1
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
    ax.set_xlabel(latexstring(xlabel), fontsize=10)
    ax.set_ylabel(latexstring(ylab), rotation=0, fontsize=10, labelpad=10)

    # Sonic line:
    if var == "M" || var == "MRel"
        ax.plot([minimum(x), maximum(x)], [1.0, 1.0], "--k")
    end

    # Marks are tagged so that FitY_ShockFlame! ignores them.
    if MarkSF
        isfinite(snap.xf) && ax.axvline(snap.xf, color="r", linestyle="--", gid="sfmark")
        isfinite(snap.xs) && ax.axvline(snap.xs, color="g", linestyle="--", gid="sfmark")
    end

    return v_plot

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
function SetTitle_ShockFlame!(fig, snap::SFSnapshot)
    fig.suptitle(latexstring("\\mathcal{Q}=", snap.q, ", M_i=", snap.ML,
                             ", M_t=", sprintf1("%.2f", snap.Mt),
                             ", t=", sprintf1("%.2e", snap.t)), fontsize=10)
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
    logVars::Vector{String}=String[],
    xlim=nothing, YLims::Dict=Dict{String,Tuple{Float64,Float64}}(), FitY::Bool=false,
    title::Bool=true, overlay=nothing)

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
        PlotVar_ShockFlame!(ax, snap, var, MarkSF=MarkSF, markersize=markersize; PlotVarKW...)
        var in logVars && ax.set_yscale("log")
    end

    overlay === nothing || overlay(axes)

    xlim === nothing || axes[1].set_xlim(xlim[1], xlim[2])

    for (ii, var) in enumerate(PlotVars)
        if haskey(YLims, var)
            axes[ii].set_ylim(YLims[var][1], YLims[var][2])
        elseif FitY
            FitY_ShockFlame!(axes[ii])
        end
    end

    title && SetTitle_ShockFlame!(fig, snap)

    return axes

end


#-------------------------------------------------------------------------------
# Plot functions
#-------------------------------------------------------------------------------

"""
    x1Plot_ShockFlame(SC, nb; kwargs...)

Plot `PlotVars` of result `nb` of study case `SC` (`nb=1000` is the latest
saved result), one per subplot. This is the single user entry point replacing
`x1Plot_ShockFlame`, `x1Plot_ShockFlame_v2` and
`Compare_ShockFlame_FXP_Overdriven`.

Keyword arguments:
- `Interactive`: browse the saved results. Keys `6` / `4` go to the next /
  previous result (the arrow keys are reserved by matplotlib for the toolbar
  views) and `keyLog` (default `l`) toggles the y scale of the subplot under the
  mouse between linear and logarithmic. The x zoom is kept when changing
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
  compared with the ZND `uRel`;
  variables without ZND counterpart are plotted without it (with a warning).
- `RTcross`: draw the crossover temperature `RT_cross` (same formula as
  `compute_FXP_Tcross`, for the `M_t` of each result) as a red dashed horizontal
  line in the `RT` subplots; key `keyRTcross` (default `T`, also `t`) switches it
  on and off. It is computed for any `M_t` (with a warning if `M_t < M_CJ`).
- `MCJline`: draw `M_CJ` as a red dashed horizontal line in the Mach number
  subplots (`MRel`); key `keyMCJ` (default `M`, also `m`) switches it on and off.
  `xf_znd`, `Deltax_znd`: length and step of the ZND solution.
- `logVars`: variables drawn with a logarithmic y axis (e.g. `["Y","p"]`).
- `FitY`: fit the y axis to the visible data. Defaults to `Interactive`.
- `YLims`: fixed y limits per variable, e.g. `Dict("Y"=>(1e-12, 2.0))`.
- `MarkSF`, `markersize`, `mFig`, `nFig`, `w`, `h`, `Deltat`: as in the
  original functions.
- `SaveFig`: save `\$(FigUbi)SC<SC>_<nb>_x1Plot.png` (`..._x1Plot_ZND.png` with
  `ZND=true`) for every drawn result.

See `PlotVar_ShockFlame!` and `RenderSnapshot_ShockFlame!` for details.
Returns `nothing`.
"""
function x1Plot_ShockFlame(SC::Int, nb::Int;
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
    RTcross::Bool=false, keyRTcross::String="T", MCJline::Bool=false, keyMCJ::String="M")

    ZND && CheckZNDVars_ShockFlame(PlotVars)
    zndOn        = Ref(ZND)
    rtcOn        = Ref(RTcross)
    mcjOn        = Ref(MCJline)
    fit          = FitY === nothing ? Interactive : FitY
    logSet       = Set(logVars)
    nb_ref       = Ref{Int}(ResolveNb_ShockFlame(SC, nb))
    first_render = Ref(true)
    busy         = Ref(false)
    home_xlims   = Ref(Any[])
    fig          = NewFigure_ShockFlame(mFig, nFig, w=w, h=h)

    function render_plot(current_nb::Int)

        if !isfile(GetFileName(SC, current_nb))
            @warn "File does not exist: $(GetFileName(SC, current_nb))"
            return false
        end

        # x view to keep (only if a previous render exists).
        prev_xlims  = first_render[] ? nothing : [ax.get_xlim() for ax in fig.axes]

        snap        = LoadSnapshot_ShockFlame(SC, current_nb, Deltat=Deltat)

        # The ZND solution is only computed and drawn for overdriven shocks (M_t > M_CJ).
        useZND      = zndOn[] && snap.Mt > MCJ(snap.GasModel.gamma, snap.q)
        if zndOn[] && !useZND
            @info "ZND not drawn for nb=$(current_nb): M_t=$(round(snap.Mt, digits=3)) <= M_CJ"
        end
        znd         = useZND ? ZND_ShockFlame(snap, xf_znd=xf_znd, Deltax_znd=Deltax_znd) : nothing
        dynKW       = useZND ? (color="r", legend=false, labelsuffix="~(\\mathrm{dyn})",
                                curvelabel="dyn", xlabel="x") : NamedTuple()
        rtc         = rtcOn[] ? RTCross_ShockFlame(snap) : nothing
        overlay     = (useZND || rtc !== nothing || mcjOn[]) ? (axs -> begin
                        useZND && OverlayZND_ShockFlame!(axs, PlotVars, znd)
                        rtc === nothing || OverlayRTCross_ShockFlame!(axs, PlotVars, rtc)
                        mcjOn[] && OverlayMCJ_ShockFlame!(axs, PlotVars, MCJ(snap.GasModel.gamma, snap.q))
                      end) : nothing
        axes        = RenderSnapshot_ShockFlame!(fig, snap, PlotVars, mFig=mFig, nFig=nFig,
                        MarkSF=MarkSF, markersize=markersize, logVars=collect(logSet),
                        YLims=YLims, FitY=fit, PlotVarKW=dynKW,
                        overlay=overlay)

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
            for (idx, ax) in enumerate(axes)
                haskey(YLims, PlotVars[idx]) || (fit && FitY_ShockFlame!(ax))
            end
            toolbar !== nothing && toolbar.push_current()
        end

        if SaveFig
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
        if toggle || toggleZ || toggleT || toggleM || (delta != 0 && target >= 1)
            busy[] = true
            try
                if toggleM
                    mcjOn[] = !mcjOn[]
                    render_plot(nb_ref[])
                elseif toggleT
                    rtcOn[] = !rtcOn[]
                    render_plot(nb_ref[])
                elseif toggleZ
                    zndOn[] = !zndOn[]
                    zndOn[] && CheckZNDVars_ShockFlame(PlotVars)
                    render_plot(nb_ref[])
                elseif toggle
                    # Toggle the y scale of the subplot under the mouse.
                    inax    = event.inaxes
                    idx     = inax === nothing ? nothing : findfirst(a -> a == inax, fig.axes)
                    if idx !== nothing && idx <= length(PlotVars)
                        var = PlotVars[idx]
                        var in logSet ? delete!(logSet, var) : push!(logSet, var)
                        render_plot(nb_ref[])
                    end
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

    fig.canvas.mpl_connect("key_press_event", on_key)

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

    xv_znd      = xv_znd .+ snap.xs
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
function OverlayZND_ShockFlame!(axes, PlotVars::Vector{String}, znd)

    c   = znd.curves
    for (ii, var) in enumerate(PlotVars)
        ax = axes[ii]
        if var == "Y"
            ax.plot(znd.x, c["YF"], "--b", gid="znd", label=latexstring("Y_F~(\\mathrm{znd})"))
            ax.plot(znd.x, c["YX"], "--g", gid="znd", label=latexstring("Y_X~(\\mathrm{znd})"))
            ax.plot(znd.x, c["YP"], "--r", gid="znd", label=latexstring("Y_P~(\\mathrm{znd})"))
        elseif var == "mdot"
            ax.plot(znd.x, c["mdotF"], "--r", gid="znd", label=latexstring("-\\dot{m}_F~(\\mathrm{znd})"))
            ax.plot(znd.x, c["mdotX"], "--g", gid="znd", label=latexstring("\\dot{m}_X~(\\mathrm{znd})"))
            ax.plot(znd.x, c["mdotP"], "--b", gid="znd", label=latexstring("\\dot{m}_P~(\\mathrm{znd})"))
        else
            name = ZNDName_ShockFlame(var)
            name === nothing && continue
            # Same color as the simulation curve of this subplot (red, see dynKW).
            ax.plot(znd.x, c[name], "--r", gid="znd", label="znd")
        end
        Legend_ShockFlame!(ax)
    end

end

"""
    Video_ShockFlame(SC; kwargs...)

Modular version of `Video_ShockFlame`: plays results `n1:min(n2, latest)` of
study case `SC`. Space or `p` pauses/resumes, `r` restarts; when the video ends
the last frame stays until `r` is pressed or the window is closed.

Extra keyword arguments with respect to the original: `logVars` (default
empty: all subplots start in linear scale), `YLims` (default `Y` in `[1e-12, 2]`,
applied only while that subplot is in log scale) and `FitY` (fit the y axis
of the other subplots to the visible data, default `true`). `TrackSF` keeps
the x axis between the shock and the flame. Key `keyLog` (default `l`) toggles
the y scale of the subplot under the mouse between linear and logarithmic; the
choice is kept in the following frames. `RTcross` (default `false`) draws `RT_cross`
in the `RT` subplots; key `keyRTcross` (default `T`) switches it on and off.
`MCJline` (default `false`) draws `M_CJ` in the `MRel` subplots; key `keyMCJ`
(default `M`) switches it on and off.
"""
function Video_ShockFlame(SC::Int;
    PlotVars::Vector{String}=["Y", "p", "vxRel", "mdot", "RT", "MRel"],
    MarkSF::Bool=false, TrackSF::Bool=true,
    mFig::Int=2, nFig::Int=3,
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50,
    tPlot::Float64=30.0, n1::Int=0, n2::Int=1000, Deltat::Float64=1.0,
    logVars::Vector{String}=String[], keyLog::String="l",
    RTcross::Bool=false, keyRTcross::String="T", MCJline::Bool=false, keyMCJ::String="M",
    YLims::Dict=Dict("Y"=>(1e-12, 2.0)), FitY::Bool=true, markersize::Real=0.2)

    nbMax       = min(load("$(ResUbi)LIRKHyp_SC$(SC)_info.jld2", "nb_SaveRes"), n2)
    history     = LoadHistory_ShockFlame(SC, nbMax, Deltat=Deltat)
    tPlot_nFig  = tPlot/(nbMax-n1)

    fig         = NewFigure_ShockFlame(mFig, nFig, w=w, h=h)
    axes        = fig.axes[1:length(PlotVars)]
    ShareX_ShockFlame!(axes)

    is_paused    = Ref(false)
    restart_flag = Ref(false)
    logSet       = Set(logVars)
    last_snap    = Ref{Any}(nothing)
    rtcOn        = Ref(RTcross)
    mcjOn        = Ref(MCJline)

    function draw_frame(snap)
        rtc     = rtcOn[] ? RTCross_ShockFlame(snap) : nothing
        RenderSnapshot_ShockFlame!(fig, snap, PlotVars, axes=axes,
            MarkSF=MarkSF, markersize=markersize, logVars=collect(logSet),
            xlim=TrackSF ? (snap.xs-0.2, snap.xf+0.2) : nothing,
            YLims=filter(kv -> kv.first in logSet, YLims), FitY=FitY,
            overlay=(rtc === nothing && !mcjOn[]) ? nothing : (axs -> begin
                rtc === nothing || OverlayRTCross_ShockFlame!(axs, PlotVars, rtc)
                mcjOn[] && OverlayMCJ_ShockFlame!(axs, PlotVars, MCJ(snap.GasModel.gamma, snap.q))
            end))
    end

    fig.canvas.mpl_connect("key_press_event", event -> begin
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
        elseif event.key == "r"
            restart_flag[] = true
            is_paused[] = false
        elseif event.key == keyLog
            # Toggle the y scale of the subplot under the mouse; it is kept in
            # the following frames, and the current frame is redrawn.
            inax    = event.inaxes
            idx     = inax === nothing ? nothing : findfirst(a -> a == inax, axes)
            if idx !== nothing
                var = PlotVars[idx]
                var in logSet ? delete!(logSet, var) : push!(logSet, var)
                if last_snap[] !== nothing
                    draw_frame(last_snap[])
                    fig.canvas.draw_idle()
                end
            end
        end
    end)

    kSave       = 0
    while true
        restart_flag[] = false

        for nb=n1:nbMax

            t0      = time()
            snap    = LoadSnapshot_ShockFlame(SC, nb, history=history)

            if !(isfinite(snap.xs) && isfinite(snap.xf) && isfinite(snap.xsdot))
                @warn "Video_ShockFlame: nb=$nb skipped (xs=$(snap.xs), xf=$(snap.xf), xsdot=$(snap.xsdot))"
                continue
            end

            # Exit if the window has been closed.
            PyPlot.plt.fignum_exists(fig.number) || return

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

        # Keep the last frame until 'r' is pressed or the window is closed.
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


function Track_ShockFlame(SC::Int; SaveFig::Bool=false, w::Float64=10.0, h::Float64=10.0, 
    Deltat::Float64=1.0, nbMax::Int=-1)

    if nbMax<0
        nbMax   = load("$(ResUbi)LIRKHyp_SC$(SC)_info.jld2", "nb_SaveRes")
    end
    
    FileName    = GetFileName(SC, nbMax)
    q           = load(FileName, "q")
    ML          = load(FileName, "ML")
    tv          = load(FileName, "tv")
    xsv         = load(FileName, "xsv")
    xfv         = load(FileName, "xfv")
    gamma       = load(FileName, "gamma")
#     Mtv         = load(FileName, "Mrelv")

    #Relative Mach number:
    u0              = ML*sqrt(gamma*(1.0+q))    #u0=uL, RT0=1, RTL=RT0*(1+q)
    a0              = sqrt(gamma*1.0)
    xsdotv          = xdotfun(tv, xsv, Deltat=Deltat)
    Mtv             = @. (u0-xsdotv)/a0
#     Mtv[xfv.<xsv]   .= NaN
    
    #Plot xs(t) and xf(t):
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.0, left=1.5, right=0.2)
    plot(tv, xsv, "b")
    plot(tv, xfv, "r")
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring("t"), fontsize=10)
    ylabel(latexstring("x_1"), fontsize=10, rotation=0)    
    legend(["shock", "flame"], fontsize=8)
#     title(latexstring("\\mathcal{Q}=", q, ", \\mathcal{M}_i=", ML), fontsize=10)
    title(latexstring("\\mathcal{Q}=", q, ", M_L=", ML), fontsize=10)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_TrackShockFlame.png", dpi=800, pad_inches=0)
    end
    
    #Plot xsdot:
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.0, left=1.5, right=0.2)
    xfdotv          = xdotfun(tv, xfv, Deltat=Deltat)
    plot(tv, xsdotv, "b")
    plot(tv, xfdotv, "r")
    legend([latexstring("\\dot{x}_s"), latexstring("\\dot{x}_f")], fontsize=8)
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring("t"), fontsize=10)
    #         ylabel(latexstring("x_f-x_s"), fontsize=10, rotation=0)    
    ylabel(latexstring("\\dot{x}"), fontsize=10, rotation=0)
    allv            = vcat(xsdotv, xfdotv)
    ylim(-0.1, 1.2*maximum(allv[.!isnan.(allv)]))
    title(latexstring("\\mathcal{Q}=", q, ", \\mathcal{M}_i=", ML), fontsize=10)
    grid("on")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_TrackShockFlame2.png", dpi=800, pad_inches=0)
    end
    
    #Plot xs(t) and xf(t):
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.0, left=1.5, right=0.2)
    auxv            = xfv-xsv
    auxv[xfv.<xsv]  .= NaN
    plot(tv, auxv, "m")
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring("t"), fontsize=10)
    ylabel(latexstring("x_f-x_s"), fontsize=10, rotation=0)    
    title(latexstring("\\mathcal{Q}=", q, ", \\mathcal{M}_i=", ML), fontsize=10)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_TrackShockFlame3.png", dpi=800, pad_inches=0)
    end
    
    #Plot Mrel(t):
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.0, left=1.5, right=0.2)
    plot(tv, Mtv, "b")
    plot([tv[1], tv[end]], fill(MCJ(gamma, q), 2), "--r")
    legend([latexstring("\\mathcal{M}_t"), latexstring("\\mathcal{M}_{CJ}")], fontsize=8)
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring("t"), fontsize=10)
    ylabel(latexstring("\\mathcal{M}"), fontsize=10, rotation=0, labelpad=10.0)    
    title(latexstring("\\mathcal{Q}=", q, ", \\mathcal{M}_i=", ML), fontsize=10)
    ylim(-0.1, 1.2*max(MCJ(gamma, q), maximum(Mtv[.!isnan.(Mtv)])))
    grid("on")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Mt.png", dpi=800, pad_inches=0)
    end
#     println("M_rel=", Mtv)
    
    return
    
end






function MtPlot(q::Float64, beta::Float64; 
    SaveFig::Bool=false, w::Float64=10.0, h::Float64=10.0)

    Miv     = []
    Mt1v    = []
    Mt2v    = []
    gamma   = 1.4
    if q==5.0 && beta==10.0
        Miv     = [ 2.0,    2.4,    2.8,    2.86,   2.92, 
                    2.99,   3.05,   3.11,   3.18,   3.26, 
                    3.33,   3.4, 
                    4.0,    4.9,    5.8,    6.9,    8.4, 
                    10.0    ]
        Mt1v    = [ 2.58,   3.26,   3.95,   4.05,   4.16, 
                    4.29,   4.41,   5.55,   5.60,   5.69, 
                    5.78,   5.86, 
                    6.65,   7.93,   9.27,   10.94,  13.24, 
                    15.71   ] 
        Mt2v    = copy(Mt1v)
    else
        error("Undefined case")
    end
    #OLD CASES:
    #=
    if q==5.0 && beta==25.0
        Miv     = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 10.0]
        Mt1v    = [2.6, 5.2, 5.4, 5.99, 6.65, 8.08, 15.71 ]
        Mt2v    = [2.6, 5.2, 5.4, 5.99, 6.65, 8.08, 15.72 ]
        #Case M1=2.5 is not totally stationary
    elseif q==5.0 && beta==50.0
        Miv     = [2.0, 2.5, 3.0, 3.5, 4.0]
        Mt1v    = [2.59, 3.42, NaN, 5.7, 5.5 ]
        Mt2v    = [2.59, 3.42, NaN, 6.3, 8.5 ]
    elseif q==1.0 && beta==25.0
        Miv     = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 10.0, 20.0]
        Mt1v    = [NaN, 3.17, 3.73, 4.30, 4.88, 6.07, 12.09, 24.13 ]
        Mt2v    = [NaN, 3.22, 3.76, 4.33, 4.91, 6.09, 12.10, 24.16 ]
    else
        error("Undefined case")
    end
    =#
    
    #Plot xs(t) and xf(t):
    PyPlotFigure(w=w, h=h, top=1.0, bottom=1.5, left=2.2, right=0.2)
    semilogx(Miv, Mt1v./Miv, "sb", markersize=3.0, label="_")
    semilogx(Miv, Mt2v./Miv, "sb", markersize=3.0, label="_nolegend_")
    Miv_ex          = logspace(minimum(Miv), maximum(Miv), 100)
    Mtv_ex          = zeros(size(Miv_ex))
    for ii=1:length(Miv_ex)
        Mtv_ex[ii], = MtMrDeton(Miv_ex[ii], gamma, q) 
    end
    semilogx(Miv_ex, Mtv_ex./Miv_ex, color=:orange, label="detonation")
    for ii=1:length(Miv_ex)
        Mtv_ex[ii], = MtMrShock(Miv_ex[ii], gamma, q) 
    end
    semilogx(Miv_ex, Mtv_ex./Miv_ex, color=:cyan, label="shock")
    #
    tick_params(axis="both", which="both", labelsize=TickSize)
#     xlabel(latexstring("\\mathcal{M}_i"), fontsize=10)
#     ylabel(latexstring("\\mathcal{M}_t/\\mathcal{M}_i"), fontsize=10, rotation=0, labelpad=20.0)
#     title(latexstring("\\mathcal{Q}=", q, ", T_I/T_s^{CJ}=T_B/T_s^{CJ}=", beta), fontsize=10)
    xlabel(latexstring("M_L"), fontsize=10)
    ylabel(latexstring("M_T/M_L"), fontsize=10, rotation=0, labelpad=20.0)    
    title(latexstring("\\mathcal{Q}=", q), fontsize=10)
    legend(fontsize=8)
    if SaveFig
        savefig("$(FigUbi)Mt_q$(q)_beta$(beta).png", dpi=800, pad_inches=0)
    end
    
end

function MCJ_(gamma::Float64, q::Float64)

    return sqrt((gamma+1)*q/2)+sqrt(1.0+(gamma+1)*q/2)
    
end

#Based on Napierialski et al (AMM2024):
function MtMrDeton(Mi::Float64, gamma::Float64, q::Float64)

    Rfun(M)     = ((gamma+1)*M^2)/((gamma-1)*M^2+2);
    Pfun(M)     = (2*gamma*M^2-gamma+1)/(gamma+1);
    
    #qq=Q, q=QQ, in Cesar's code
#     qq          = 0.0       
    qq          = q
    
    #Define residual function:
    Ri          = Rfun(Mi)
    Ai          = sqrt(Pfun(Mi)/Ri)
    function ResidualFun!(x::Vector{Float64}, f::Union{Vector{Float64},Nothing}, g::Union{Vector{Float64},Nothing})
        
#         println("Mt=",x[1])
        
        #Transmitted shock:
        Mt      = x[1]
        if Mt<0
            @warn("Negative Mach")
            return -1
        end
        kappa2  = ((Mt^2-1)^2 - 2*(gamma+1)*qq*Mt^2)/Mt^4
        if kappa2<0
            @warn("Negative kappa2")
            return -1
        end
        kappa   = sqrt(kappa2)
        Rt      = (gamma+1)*Mt^2 / ((gamma-kappa)*Mt^2 +1)
        
        #Reflected shock:
        Mr      = sqrt( (3 - gamma - 2*Mi^2 + 2*gamma*Mi^2 + (1+gamma)*(1+kappa)*Mt^2) / 
                        (2*(1+gamma*(2*Mi^2-1))) )
        Rr      = Rfun(Mr)
      
        #Residual:
        f[1]    = Mi*(1.0-1.0/Ri)-Ai*Mr*(1.0-1.0/Rr)-Mt*(1.0-1.0/Rt)/sqrt(1.0+q)    #Eq. (6)
        if !isnothing(g)
            faux    = zeros(1)
            ResidualFun!([Mt+1e-4], faux, nothing)
            J       = (faux[1]-f[1])/1e-4
            J       *= 10.0      #Safety factor to trick Anderson's method
            g[1]    = J\f[1]
        end
        
        return 1
        
    end
    
    #DEBUG:
    if false
        Nplot       = 100
        fplotv      = zeros(Nplot)
        gplotv      = zeros(Nplot)
        Mtplotv     = linspace(0.9*Mi, MCJ_(gamma, q)+Mi, Nplot);
        for ii=1:length(Mtplotv)
            faux    = zeros(1)
            gaux    = zeros(1)
            flag    = ResidualFun!([Mtplotv[ii]], faux, gaux)
            if flag<0
                fplotv[ii]  = NaN
                gplotv[ii]  = NaN
            else
                fplotv[ii]  = faux[1]
                gplotv[ii]  = gaux[1]
            end
        end
        figure()
        plot(Mtplotv, fplotv, "b")
        plot(Mtplotv, gplotv, "r")
        axis("equal")
        grid("on")
        display(MCJ_(gamma, q) + Mi)
        println("")
    #     return
    end
    
    #Initial conditions:
    Mcj         = MCJ_(gamma, q)
    Mt0         = Mi+Mcj
    
    #Call Anderson's method:
    NLS_output  = Anderson(FW_NLS((u,g)->ResidualFun!(u,[0.0],g)), 
                    [Mt0], 
                    AbsTolX=0.0, RelTolX=0.0, 
                    AbsTolG=1e-4, RelTolG=0.0,
                    Display="notify", MaxIter=100, memory=0)
    Mt          = NLS_output[1][1]
    kappa       = sqrt(((Mt^2-1)^2 - 2*(gamma+1)*qq*Mt^2)/Mt^4)
    Mr          = sqrt( (3 - gamma - 2*Mi^2 + 2*gamma*Mi^2 + (1+gamma)*(1+kappa)*Mt^2) / 
                        (2*(1+gamma*(2*Mi^2-1))) )
    if NLS_output[2].flag<0
        Mt      = NaN
        Mr      = NaN
    end
    
    return Mt, Mr
    
end

#Based on Napierialski et al (AMM2024):
function MtMrShock(Mi::Float64, gamma::Float64, q::Float64)

    Rfun(M)     = ((gamma+1)*M^2)/((gamma-1)*M^2+2);
    Pfun(M)     = (2*gamma*M^2-gamma+1)/(gamma+1);
    
    #Define residual function:
    Ri          = Rfun(Mi)
    Ai          = sqrt(Pfun(Mi)/Ri)
    function ResidualFun!(x::Vector{Float64}, f::Union{Vector{Float64},Nothing}, g::Union{Vector{Float64},Nothing})
    
        Mt      = x[1]
        Mr      = x[2]
        Rt      = Rfun(Mt)
        Rr      = Rfun(Mr)
        f[1]    = Mi*(1.0-1.0/Ri)-Ai*Mr*(1.0-1.0/Rr)-Mt*(1.0-1.0/Rt)/sqrt(1.0+q)            #Eq. (6)
        f[2]    = Mr^2*(1.0+gamma*(2.0*Mi^2-1)) - Mt^2*(gamma+1)-Mi^2*(gamma-1)+(gamma-1)   #Eq. (7)
        if !isnothing(g)
            J       = [ -(1-1/Rt)/sqrt(1+q)     -Ai*(1-1/Rr)    
                        -2*Mt*(gamma+1)         2*Mr*(1+gamma*(2*Mi^2-1)) ]
            g       .= J\f
#             g       .= f
        end
        
        return 1
        
    end
    
#     fv      = zeros(2)
#     ResidualFun!([2.58,1.26], fv, nothing)
#     display(fv)
#     return
    
    #Initial conditions:
    Mt0         = Mi
    Mr0         = sqrt( (Mt0^2*(gamma+1)+Mi^2*(gamma-1)+(gamma-1))/(1.0+gamma*(2.0*Mi^2-1)) )
    
    #Call Anderson's method:
    NLS_output  = Anderson(FW_NLS((u,g)->ResidualFun!(u,[0.0, 0.0],g)), 
                    [Mt0, Mr0], 
                    AbsTolX=0.0, RelTolX=0.0, 
                    AbsTolG=1e-4, RelTolG=0.0,
                    Display="notify", MaxIter=100)
    Mt          = NLS_output[1][1]
    Mr          = NLS_output[1][2]
    if NLS_output[2].flag<0
        Mt      = NaN
        Mr      = NaN
    end
    
    return Mt, Mr
    
end

#Compute shock and flame velocities:
function xdotfun(tv::Vector{Float64}, xv::Vector{Float64};
    Deltat::Float64=1.0)
#     return xdotfun_nofilter(tv, xv)
#     return xdotfun_filter1(tv, xv)
#     return xdotfun_filter2(tv, xv)
#     return xdotfun_filter3(tv, xv)
    return xdotfun_filter4(tv, xv; Deltat=Deltat)
end

function xdotfun_nofilter(tv::Vector{Float64}, xv::Vector{Float64})

    N       = length(tv)
    if N==1
        return [ 0.0 ]
    end
    
    xdotv           = zeros(N)
    xdotv[1]        = (xv[2]-xv[1])/(tv[2]-tv[1])
    xdotv[2]        = (xv[2]-xv[1])/(tv[2]-tv[1])
    for ii=3:N
        xdotv[ii]   = (3*xv[ii]-4*xv[ii-1]+xv[ii-2])/(3*tv[ii]-4*tv[ii-1]+tv[ii-2])
#         xdotv[ii]   = (xv[ii]-xv[ii-1])/(tv[ii]-tv[ii-1])
    end
    
    return xdotv
    
end

function xdotfun_filter1(tv::Vector{Float64}, xv::Vector{Float64})

    N               = length(tv)
    if N==1
        return [ 0.0 ]
    end
    
    #Evaluate mean position at t^(n+1/2):
    t12v            = 0.5*(tv[1:N-1]+tv[2:N])
    x12v            = 0.5*(xv[1:N-1]+xv[2:N])
    
    #Velocities from linear interpolation of mean position:
    xdotv           = zeros(N)
    xdotv[1]        = NaN
    for ii=2:N-1
        xdotv[ii]   = (x12v[ii]-x12v[ii-1])/(t12v[ii]-t12v[ii-1])
    end
    xdotv[N]        = NaN
    
    return xdotv
    
end

function xdotfun_filter2(tv::Vector{Float64}, xv::Vector{Float64}; 
    DeltaN::Int=10, Nfilter::Int=1)
    
    N       = length(tv)
    if N==1
        return [ 0.0 ]
    end
    
    #Compute velocities with second order formula. The result is noisy:
    xdotv           = zeros(N)
    xdotv[1]        = (xv[2]-xv[1])/(tv[2]-tv[1])
    xdotv[2]        = (xv[2]-xv[1])/(tv[2]-tv[1])
    for ii=3:N
        xdotv[ii]   = (3*xv[ii]-4*xv[ii-1]+xv[ii-2])/(3*tv[ii]-4*tv[ii-1]+tv[ii-2])
    end
    xdot0v          = 1.0*xdotv
    
    #Filter velocities:
    for ifilter=1:Nfilter
        mass            = 0.0
        inertia1        = 0.0
        tbar            = 0.0
        xdothatv        = 0.0*xdotv
        for ii=2:1:N
        
            #Update length and barycenter and inertia:
            i0              = max(1, ii+1-DeltaN)
            L_new           = tv[ii]-tv[i0]
            tbar_new        = 0.5*(tv[i0]+tv[ii])
            inertia1_new    = inertia1 + mass*(tbar-tbar_new)
            mass_new        = mass
            
            #Update mass and momentum respect to center of time interval [t_(i-1), t_i] 
            #with new element:
            h               = tv[ii]-tv[ii-1]
            Delta_mass      = 0.5*h*xdotv[ii-1] + 0.5*h*xdotv[ii]
            mass_new        += Delta_mass
            inertia1_new    += h*h/12*(xdotv[ii]-xdotv[ii-1]) + 
                                Delta_mass*(0.5*(tv[ii]+tv[ii-1])-tbar_new)
            
            #Update mass and momentum respect to center of time interval [t_(i-1), t_i] 
            #with outgoing element:
            if ii>DeltaN
                h               = tv[ii-DeltaN+1]-tv[ii-DeltaN]
                Delta_mass      = 0.5*h*xdotv[ii-DeltaN+1] + 0.5*h*xdotv[ii-DeltaN]
                mass_new        -= Delta_mass
                inertia1_new    -= h*h/12*(xdotv[ii-DeltaN+1]-xdotv[ii-DeltaN]) + 
                                    Delta_mass*(0.5*(tv[ii-DeltaN+1]+tv[ii-DeltaN])-tbar_new)
            end
            
            #Seek linear xdothat such that mass and momentum are the same in
            #the interval [tv[1], tv[ii]]:
            #   0.5*L*v1 + 0.5*L*v2     = mass_new
            #   -L^2/12*v1 + L^2/12*v2  = inertia1_new
            #v1 and v2 are the values of xdothat at tv[1] and tv[ii]
            #   v1 + v2                 = mass_new/(0.5*L)
            #   -v1 + v2                = inertia1_new*12/L^2
            xdothatv[ii]    = 0.5*( mass_new/(0.5*L_new) + inertia1_new*12/L_new^2 )
        
            #Update mass and inertia1:
            mass            = mass_new
            inertia1        = inertia1_new
            tbar            = tbar_new
            
        end
        xdotv               .= xdothatv
    end
    
    #NaNs for first time steps:
    xdotv[1:DeltaN]     .= NaN
        
#     figure()
#     plot(tv, xdot0v, "x-g")
#     plot(tv, xdotv, "+-r")
#     ylim(-0.5, 2.5)
#     error("")
    
    return xdotv
    
end

function xdotfun_filter3(tv::Vector{Float64}, xv::Vector{Float64};
    DeltaN::Int=50)
    
    #Compute velocities with first order formula. Compute average of 
    #DeltaN last time steps:
    N               = length(xv)
    xdotv           = zeros(N)
    xdotv[1]        = NaN
    for ii=2:N
        jmax            = min(DeltaN,ii-1)
        for jj=1:jmax
            xdotv[ii]   += (xv[ii]-xv[ii-jj])/(tv[ii]-tv[ii-jj])
        end
        xdotv[ii]       /= jmax
    end
        
    return xdotv
    
end

function xdotfun_filter4(tv::Vector{Float64}, xv::Vector{Float64};
    Deltat::Float64=1.0)
    
    #Compute velocities with first order formula. Take a minimum time step
    #given by Deltat:
    N               = length(xv)
    xdotv           = zeros(N)
    xdotv[1]        = NaN
    for ii=2:N
        
        #Find maximum jj such that tv[ii]-tv[jj]>Deltat
        jj          = ii-1
        while true
            if tv[ii]-tv[jj]>Deltat
                break
            end
            if jj==1
                break
            end
            jj      -= 1
        end
            
        #Compute velocity:
        xdotv[ii]   = (xv[ii]-xv[jj])/(tv[ii]-tv[jj])
        
    end
        
    return xdotv
    
end

#NOTAS:
#   -¿Se podría definir de forma más rigurosa la posición de la llama y de la onda?
#   -¿Sería necesario hacerlo para reducir las oscilaciones?
#   -Faltaría modificar el código para que la onda esté quieta.
#   -Con beta=50.0 los resultados no parecen ser estacionarios.
#   -Hay casos en los que la llama se separa de la onda al principio y luego
#   se junta (utilizando el modelo de Arrhenius). Sería importante limitar el dominio
#   o el tiempo.

#=
function pYref_ShockFlame(SC::Int)

    solver      = GetSolver(SC,0)
    GasModel    = GasFXP_ShockFlame(SC,0)
    
    pv          = SolutionAtNodes(solver,GasModel,"p")
    YFv         = SolutionAtNodes(solver,GasModel,"Y_F")
    
    pref        = 0.5*(maximum(pv)+minimum(pv))
    YFref       = 0.5*(maximum(YFv)+minimum(YFv))
   
    return pref, YFref
    
end

function xsxf_ShockFlame(solver::SolverDataSave, GasModel::GasFXP,
    pref::Float64, YFref::Float64)

    #Pressure and YF at mesh nodes:
    x1v         = solver.mesh.NodesCoords[:,1]
    pv          = SolutionAtNodes(solver,GasModel,"p")
    YFv         = SolutionAtNodes(solver,GasModel,"Y_F")
        
    #Sort by x1 coordinate:
    sortv       = sortperm(x1v)
    x1v         = x1v[sortv]
    pv          = pv[sortv]
    YFv         = YFv[sortv]
    
    #Get shock and flame positions:
    xf          = NaN
    T           = findlast(YFv.>YFref)
    if T<length(x1v)
        xf      = x1v[T] + (x1v[T+1]-x1v[T])*(YFref-YFv[T])/(YFv[T+1]-YFv[T])
    end
    #
    xs          = NaN
    T           = findlast(pv.<pref)
    if T<length(x1v)
        xs      = x1v[T] + (x1v[T+1]-x1v[T])*(pref-pv[T])/(pv[T+1]-pv[T])
    end
        
    return xs, xf
    
end

=#

function Table_ShockFlame(SC::Int, nb::Int;
    PlotVars::Vector{String}=["x", "Y_F", "Y_X", "Y_P", "vxRel", "p", "RT", "rho", "MRel"],
    MarkSF::Bool=false,
    mFig::Int=3, nFig::Int=2, 
    SaveFig::Bool=false, w::Float64=9.50, h::Float64=6.50, 
    Deltat::Float64=1.0)
    
    #Extract data:
    if nb==1000
        nb      = load("$(ResUbi)LIRKHyp_SC$(SC)_info.jld2", "nb_SaveRes")
    end
    GasModel    = GasFXP_ShockFlame(SC,nb)
    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    q           = load(FileName, "q")
    ML          = load(FileName, "ML")
    
    #Get positions of shock and flame
    tv          = load(FileName, "tv")
    xsv         = load(FileName, "xsv")
    xfv         = load(FileName, "xfv")
    xsdotv      = xdotfun(tv, xsv, Deltat=Deltat)
    Nt          = length(tv)
    xs          = xsv[Nt]
    xf          = xfv[Nt]
    xsdot       = xsdotv[Nt]
    
    #Create matrix:
    nNodes      = solver.fes.PSpace.nDof
    A           = zeros(nNodes, length(PlotVars))
    
    #No zoom - loop plot variables:
    for ii=1:length(PlotVars)
        
        if PlotVars[ii]=="vxRel" #unfinished - xsxdot not defined
            vxRel           = SolutionAtNodes_vxRel(solver, GasModel, xsdot)
            A[:,ii]         = vxRel
        elseif PlotVars[ii]=="MRel"
            MRel            = SolutionAtNodes_MRel(solver, GasModel, xsdot)
            A[:,ii]         = MRel
        elseif PlotVars[ii]=="x"
            A[:,ii]         = solver.fes.PSpace.NodesCoords[:,1]
        else    
            A[:,ii]         = PlotNodes((x1,x2)->x1, solver, GasModel, PlotVars[ii])
        end
        
    end
    
    #Sort by x-coord:
    xv          = solver.fes.PSpace.NodesCoords[:,1]
    permv       = sortperm(xv)
    A           = A[permv,:]    
        
    if SaveFig
        open("$(FigUbi)SC$(SC)_$(nb).csv"; write=true) do f
            for ii=1:length(PlotVars)-1
                write(f, PlotVars[ii], '\t')
            end
            for ii=length(PlotVars)
                write(f, PlotVars[ii], '\n')
            end
            writedlm(f, A, '\t')
        end
#         writedlm(, A, '\t')
    end
        
    return
    
end