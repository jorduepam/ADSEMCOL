include("../../../src/AuxiliaryFunctions/ODESolvers.jl")

"""
    ChapmanJouguet(gamma, q) -> M0, rho0, u0, p0, rho1, u1, p1, rho2, u2, p2

Chapman-Jouguet detonation states for a perfect gas with heat release `q`
(`rho0*Q0*YF0 = q*gamma/(gamma-1)*p0`), where
  - state 0: unburnt gas ahead of the shock (normalized to rho0=p0=1),
  - state 1: post-shock (von Neumann) state,
  - state 2: burnt gas at the CJ (sonic) point.
`M0` is the CJ Mach number of the detonation.
"""
function ChapmanJouguet(gamma::Float64, q::Float64)

    #We assume the state at the end p2=1, rho2=1, u2=sqrt(gamma), YF2=0.
    p2      = 1.0
    rho2    = 1.0
    u2      = sqrt(gamma*p2/rho2)
    
    #Then, we compute the state "0" by imposing conservation equations:
    #   rho0*u0                                                 = M
    #   p0+rho0*u0^2                                            = P
    #   (gamma/(gamma-1)*p0 + 0.5*rho0*u0^2 + rho0*YF0*Q0)*u0   = HU
    #Note that rho0*Q0*YF0 = q*gamma/(gamma-1)*p0, so the last equation is
    #   ((1+q)*gamma/(gamma-1)*p0 + 0.5*rho0*u0^2)*u0           = HU
    #If we write the equations in terms of u0:
    #A u0^2 + B*u0 + HU = 0
    M       = rho2*u2
    P       = p2+rho2*u2^2
    HU      = (gamma/(gamma-1.0)*p2 + 0.5*rho2*u2^2)*u2
    #
    A       = ((1.0+q)*gamma/(gamma-1.0)-0.5)*M
    B       = -(1.0+q)*gamma/(gamma-1.0)*P
    C       = HU
    #
    u0      = (-B+sqrt(B^2-4*A*C))/(2*A)
    rho0    = M/u0
    p0      = P-M*u0
    M0      = u0/sqrt(gamma*p0/rho0)
    
    #Rescale solution with entry conditions:
    rho2    *= 1.0/rho0
    p2      *= 1.0/p0
    u2      = sqrt(gamma*p2/rho2)
    #
    rho0    = 1.0
    p0      = 1.0
    u0      = M0*sqrt(gamma)
    
    #State at region 1:
    rho1            = rho0 / ( (gamma-1.0)/(gamma+1.0) + 2.0/(gamma+1.0)/M0^2 )
    u1              = rho0*u0/rho1
    p1              = p0+rho0*u0^2-rho1*u1^2
    
    #Check results:
#     display(rho0*u0 - rho2*u2)
#     display(p0+rho0*u0^2 - (p2+rho2*u2^2))
#     display(((1+q)*gamma/(gamma-1)*p0 + 0.5*rho0*u0^2)*u0 - 
#             (gamma/(gamma-1)*p2 + 0.5*rho2*u2^2)*u2)
#     display(rho1*u1 - rho2*u2)
#     display(p1+rho1*u1^2 - (p2+rho2*u2^2))
#     display(((1+q)*gamma/(gamma-1)*p1 + 0.5*rho1*u1^2)*u1 - 
#             (gamma/(gamma-1)*p2 + 0.5*rho2*u2^2)*u2)
    
    
    return M0, rho0, u0, p0, rho1, u1, p1, rho2, u2, p2
    
end

"""
    M0_ChapmanJouguet(gamma, q)

Analytical CJ Mach number (minimum detonation Mach number) for given `gamma` and `q`.
"""
function M0_ChapmanJouguet(gamma::Float64, q::Float64)

    return sqrt(0.5*(gamma+1)*q+1) + sqrt(0.5*(gamma+1)*q)
    
end
    
"""
    FP_ChapmanJouguet(gamma, q, beta0; Deltax, xf, YF0, PlotRes) -> xv, rhov, uv, pv, YFv

Structure of a CJ detonation with one-step irreversible kinetics (fuel F -> products P),
marched in space with an explicit RK scheme from the post-shock state (state 1).
`beta0 = RT_a/RT_0` is the normalized activation temperature, `YF0` the initial fuel mass
fraction, `Deltax` the step and `xf` the final position. If `PlotRes=true`, the profiles
are plotted and compared with the CJ state 2.
"""
function FP_ChapmanJouguet(gamma::Float64, q::Float64, beta0::Float64; 
    Deltax::Float64=1e-3, xf::Float64=20.0, YF0::Float64=0.01, PlotRes::Bool=false)

    #Compute detonation states:
    M0, rho0, u0, p0, rho1, u1, p1, rho2, u2, p2 = ChapmanJouguet(gamma,q)
    
    #Compute constants:
    RT0             = p0/rho0
    RT1             = p1/rho1
    RTa             = beta0*RT0                     #Activation temperature
    beta1           = RTa/RT1                       #beta coeff with respect to CJ-1 state
    Q               = gamma/(gamma-1)*RT0/YF0*q     #Fuel formation mass enthalpy 
    BB              = u1*exp(beta1)/(1.0*beta1)     #Characteristic length is 1
                        
    #Function to compute rho, u, p, YF from vector
    #   F=[rhoY_F u, rho u, p + rho u^2, ( gamma/(gamma-1)*p + rhoY_F Q + 1/2 rho u^2 ) u ]
    function ExtractVars(Fv::Vector{Float64})
    
        #YF = rho YF u / (rho u):
        YF          = Fv[1]/Fv[2]
        
        #F2 = rho u
        #F3 = p + rho u^2 
        #F4 = ( gamma/(gamma-1)*p + 1/2 rho u^2 + rho Q YF ) u
        #   = gamma/(gamma-1)*p*u + 1/2 F2 u^2 + F1 Q
        #Combining these equations:
        #(gamma/(gamma-1)-1/2) F2 U^2 - gamma/(gamma-1) F3 U + F4-F1*Q = 0
        A       = (gamma/(gamma-1)-0.5)*Fv[2]
        B       = -gamma/(gamma-1)*Fv[3]
        C       = Fv[4]-Fv[1]*Q
        u       = (-B-sqrt(B^2-4*A*C))/(2*A)    #Smaller root: subsonic branch; the larger root is the supersonic state 0
        
        #Therefore:
        rho         = Fv[2]/u
        p           = Fv[3]-Fv[2]*u
        
        return rho, u, p, YF
        
    end
    
    #Function to march flux variables in space:
    function dF_dx(Fv::Vector{Float64})
    
        rho, u, p, YF   = ExtractVars(Fv)
        RT              = p/rho
        mFdot           = -BB*rho*YF*exp(-RTa/RT)
        return [ mFdot; 0; 0; 0 ]
        
    end
    
    #Initialize variables:
    x_n             = 0.0
    F_n             = [ rho0*YF0*u0 ;
                        rho0*u0;
                        p0 + rho0*u0^2;
                        (gamma/(gamma-1.0)*p0 + rho0*YF0*Q + 0.5*rho0*u0^2)*u0 ]
                        
#     display((rho0, u0, p0, YF0))
#     display((rho1, u1, p1, YF0))
#     display(ExtractVars(F_n))
#     return
    
    #Apply RK method:
    RK              = RK_Coefficients("Ascher3")
    NF              = length(F_n)
    fm              = zeros(NF, RK.stages)
    F_np1           = zeros(NF)
    xv              = [ -xf*0.1, x_n, x_n ]
    rhov            = [ rho0, rho0, rho1 ]
    uv              = [ u0, u0, u1 ]
    pv              = [ p0, p0, p1 ]
    YFv             = [ YF0, YF0, YF0 ]
    while x_n < xf
    
        #Compute next stage:
        fm[:,1]         = dF_dx(F_n)
        for ii=2:RK.stages
            S           = zeros(NF)
            for jj=1:ii-1
                S       += RK.AE[ii,jj]*fm[:,jj]
            end
            F_np1       .= F_n + Deltax*S
            fm[:,ii]    = dF_dx(F_np1)
        end
        
        #Update solution:
        x_n         += Deltax
        F_n         .= F_np1
        
        #Save solution:
        push!(xv, x_n)
        rho_n, u_n, p_n, YF_n   = ExtractVars(F_n)
        push!(rhov, rho_n)
        push!(uv, u_n)
        push!(pv, p_n)
        push!(YFv, YF_n)
        
    end
    
    #Compute Mach number:
    Mv      = @. uv/sqrt(gamma*pv/rhov)
    
    #Plot results:
    if PlotRes
    
        Nx      = length(xv)
        display("theoretical")
        display((rho2, u2, p2, 0.0, 1.0))
        display("numerical")
        display((rhov[Nx], uv[Nx], pv[Nx], YFv[Nx], Mv[Nx]))
        
        figure()
        plot(xv, rhov, "b")
        plot(xv, uv, "g")
        plot(xv, pv, "c")
        plot(xv, pv./rhov, "r")
        plot(xv, YFv./YF0, color="orange")
        plot(xv, Mv, "k")
        legend([latexstring("\\rho"), 
                latexstring("u"), 
                latexstring("p"),
                latexstring("T"), 
                latexstring("Y_{F}/Y_{F0}"), 
                latexstring("M") ], 
                loc="best")
        xlabel("x") #Here, x is the spatial coordinate
        grid("on")
        
    end
    
#     figure()
#     plot(xv, rhov.*uv.*YFv, "b")
#     plot(xv, rhov.*uv, "g")
#     plot(xv, pv+rhov.*uv.^2, "r")
#     plot(xv, gamma/(gamma-1)*pv.*uv + 0.5*rhov.*uv.^3 + rhov.*YFv.*Q.*uv, "k")
    
    return xv, rhov, uv, pv, YFv
    
end

"""
    FXP_ChapmanJouguet(; gamma, q, RThatI, RThatB, BIhat, BBhat, Deltax, xf, YF0, YX0, PlotRes)
        -> xv, rhov, uv, pv, YFv, YXv

Structure of a CJ detonation with three-species kinetics F -> X -> P, where F is fuel,
X an intermediate (radical) species and P the products:
  - F -> X by initiation (rate `kI`) and chain-branching (rate `kB`),
  - X -> P by recombination (rate `kR = BR`).
`RThatI`, `RThatB` are the activation temperatures of initiation and branching normalized
with the post-shock temperature RT_1; `BIhat`, `BBhat` are the normalized pre-exponential
factors. The ODE system is integrated in space with an IMEX RK scheme ("Ascher3") from the
post-shock state, using `Deltax` as step up to `xf`. If `PlotRes=true`, several figures
are produced (the last one with shared x-axis: adim. variables, mass fractions and mdot).
"""
function FXP_ChapmanJouguet(; gamma::Float64=1.4, q::Float64=5.0, 
    RThatI::Float64=10.0, RThatB::Float64=10.0, #RThat=RT_a/RT_1, #1: means post-shock condition
    BIhat::Float64=2.5e-3, BBhat::Float64=4.4e5, 
    Deltax::Float64=1e-3, xf::Float64=30.0, YF0::Float64=1.0-1e-6, YX0::Float64=1e-6, 
    PlotRes::Bool=false)

    #Compute detonation states:
    M0, rho0, u0, p0, rho1, u1, p1, rho2, u2, p2 = ChapmanJouguet(gamma,q)
    
    #Characteristic variables:
#     rhoc            = rho1
#     uc              = u1
#     pc              = rhoc*uc^2
#     tc              = 1.0
#     lc              = uc*tc
    
    #Scale pre and post shock conditions:
#     rho0            /= rhoc
#     u0              /= uc
#     p0              /= pc
#     rho1            /= rhoc
#     u1              /= uc
#     p1              /= pc
#     rho2            /= rhoc
#     u2              /= uc
#     p2              /= pc
    
    #Compute RT:
    RT0             = p0/rho0
    RT1             = p1/rho1
    
    #Compute constants:
    Q               = gamma/(gamma-1)*RT0/YF0*q     #F and X formation mass enthalpy 
    RTI             = RThatI*RT1                    #Activation temperature of initiation
    RTB             = RThatB*RT1                    #Activation temperature of chain-branching
    tc              = 1.0                           #Characteristic time
    BR              = 1.0/tc                        #t_c = B_R^{-1} = 1.0
    BI              = BIhat * BR
    BB              = BBhat * BR/rho1
                      
    Tcross = compute_FXP_Tcross(gamma=gamma, q=q, RThatB=RThatB, BBhat=BBhat, YF0=YF0)

    #Function to compute rho, u, p, YF, YX from vector 
    #   F   =   [   rhoY_F u, 
    #               rhoY_X u, 
    #               rho u, 
    #               p + rho u^2, 
    #               ( gamma/(gamma-1)*p + rho (Y_F+Y_X) Q + 1/2 rho u^2 ) u ]
    function ExtractVars(Fv::Vector{Float64})
    
        #YF = rho YF u / (rho u):
        YF          = Fv[1]/Fv[3]
        YX          = Fv[2]/Fv[3]
        YP          = 1-YF-YX
        
        #F3 = rho u
        #F4 = p + rho u^2 
        #F5 = ( gamma/(gamma-1)*p + 1/2 rho u^2 + rho Q (YF+YX) ) u
        #   = gamma/(gamma-1)*p*u + 1/2 F3 u^2 + (F1+F2) Q
        #Combining these equations:
        #(gamma/(gamma-1)-1/2) F3 U^2 - gamma/(gamma-1) F4 U + F5 - (F1+F2)*Q = 0
        A       = (gamma/(gamma-1)-0.5)*Fv[3]
        B       = -gamma/(gamma-1)*Fv[4]
        C       = Fv[5]-(Fv[1]+Fv[2])*Q
        if B^2-4*A*C<0
#             return NaN, NaN, NaN, NaN, NaN
        end
        u       = (-B-sqrt(B^2-4*A*C))/(2*A)    #Smaller root: subsonic branch (states 1 and 2); the larger root is the supersonic state 0
        
        #Therefore:
        rho         = Fv[3]/u
        p           = Fv[4]-Fv[3]*u
        
        return rho, u, p, YF, YX
        
    end
    
    #Function to march flux variables in space:
    function dF_dx(Fv::Vector{Float64}, ComputeJ::Bool)
    
        rho, u, p, YF, YX   = ExtractVars(Fv)
        RT                  = p/rho
        kI                  = BI*exp(-RTI/RT)
        kB                  = BB*exp(-RTB/RT)
        kR                  = BR
        mFdot               = - kI*rho*YF - kB*rho^2*YF*YX
        mXdot               =   kI*rho*YF + kB*rho^2*YF*YX - kR*rho*YX
        
        f                   = [ mFdot; mXdot; 0; 0; 0 ]
        J                   = zeros(0,0)
        if ComputeJ
            nVars           = length(Fv)
            J               = zeros(nVars, nVars)
            delta           = 1e-5
            for JJ=1:nVars
                Fpert       = copy(Fv)
                Fpert[JJ]   += delta
                fpert       = dF_dx(Fpert, false)[1]
                J[:,JJ]     = (fpert-f)/delta 
            end
        end
        
        return f, J
        
    end
    
    #Initialize variables:
    x_n             = 0.0
    F_n             = [ rho0*YF0*u0 ;
                        rho0*YX0*u0 ;
                        rho0*u0;
                        p0 + rho0*u0^2;
                        (gamma/(gamma-1.0)*p0 + rho0*(YF0+YX0)*Q + 0.5*rho0*u0^2)*u0 ]
                        
#     display((rho0, u0, p0, YF0, YX0))
#     display((rho1, u1, p1, YF0, YX0))
#     display(ExtractVars(F_n))
#     return
    
    #Apply RK method:
    RK              = RK_Coefficients("Ascher3")
    NF              = length(F_n)
    fm              = zeros(NF, RK.stages)
    JFm             = zeros(NF, RK.stages)
    F_np1           = zeros(NF)
    xv              = [ -xf*0.1, x_n, x_n ]
    rhov            = [ rho0, rho0, rho1 ]
    uv              = [ u0, u0, u1 ]
    pv              = [ p0, p0, p1 ]
    YFv             = [ YF0, YF0, YF0 ]
    YXv             = [ YX0, YX0, YX0 ]
        
    while x_n < xf
    
        #println("x_n=", x_n)
        
        #Compute next stage:
        fm[:,1], J      = dF_dx(F_n, true)
        JFm[:,1]        = J*F_n
        A_F             = factorize(eye(5)-Deltax*RK.AI[2,2]*J)
        for ii=2:RK.stages
            S           = zeros(NF)
            for jj=1:ii-1
                S       += RK.AE[ii,jj]*fm[:,jj] + (RK.AI[ii,jj]-RK.AE[ii,jj])*JFm[:,jj]
            end
            F_np1       .= A_F\(F_n + Deltax*S)
            fm[:,ii]    = dF_dx(F_np1, false)[1]
            JFm[:,ii]   = J*F_np1
        end
        
        #Update solution:
        x_n             += Deltax
        F_n             .= F_np1
        
        #Save solution:
        push!(xv, x_n)
        rho_n, u_n, p_n, YF_n, YX_n     = ExtractVars(F_n)
        push!(rhov, rho_n)
        push!(uv, u_n)
        push!(pv, p_n)
        push!(YFv, YF_n)
        push!(YXv, YX_n)
        
    end
    
    #Compute Mach number:
    Mv      = @. uv/sqrt(gamma*pv/rhov)
    
    #Plot each variable:
    if PlotRes
    
        Nx      = length(xv)
        println("Exact solution: rho2=$rho2, u2=$u2, p2=$p2, YF2=0.0, M2=1.0")
        println("Numerical solution: rho2=$(rhov[Nx]), u2=$(uv[Nx]), p2=$(pv[Nx]), ", 
                "YF2=$(YFv[Nx]), M2=$(Mv[Nx])")
        #=
        # Dimensional variables (not normalized):
        figure()
#         plot(xv/lc, rhov/rho1, "b")
#         plot(xv/lc, uv/u1, "g")
#         plot(xv/lc, pv/p1, "c")
#         plot(xv/lc, (pv./rhov)/RT1, "r")
        plot(xv, rhov, "b", linewidth=2.5)
        plot(xv, uv, "y", linewidth=2.5)
        plot(xv, pv, "g", linewidth=2.5)
        plot(xv, pv./rhov, "r", linewidth=2.5)
        plot(xv, Mv, "k", linewidth=2.5)
        legend([latexstring("\\rho"), 
                latexstring("u"), 
                latexstring("p"),
                latexstring("T"), 
                latexstring("M") ], 
                loc="best")
        xlabel("x") 
        grid("on")

        # Adim with respect to post-shock conditions:
        figure()
#         plot(xv/lc, rhov/rho1, "b")
#         plot(xv/lc, uv/u1, "g")
#         plot(xv/lc, pv/p1, "c")
#         plot(xv/lc, (pv./rhov)/RT1, "r")
        plot(-xv, rhov ./rho1, "b", linewidth=2.5) # rho/rho1
        plot(-xv, uv ./u1, "y", linewidth=2.5) # u/u1
        plot(-xv, pv ./p1, "g", linewidth=2.5) # p/p1
        plot(-xv, pv./rhov ./p1 .*rho1, "r", linewidth=2.5) # RT/RT1
        plot(-xv, Mv, "k", linewidth=2.5)
        legend([latexstring("\\rho/\\rho_1"), 
                latexstring("u/u_1"), 
                latexstring("p/p1"),
                latexstring("T/T1"), 
                latexstring("M") ], 
                loc="best")
        xlabel("x") 
        grid("on")

        
        figure()
        plot(-xv, YFv, color="b", linewidth=2.5)
        plot(-xv, YXv, color="orange", linewidth=2.5)
        plot(-xv, 1.0.-YFv.-YXv, color="g", linewidth=2.5)
        legend([latexstring("Y_{F}"), 
                latexstring("Y_{X}"),
                latexstring("Y_{P}")], 
                loc="best")
        xlabel("x") 
        grid("on")
        =#


        # Combined figure with shared x-axis: adim. variables, mass fractions and mdot
        RTv         = pv./rhov
        kIv         = @. BI*exp(-RTI/RTv)
        kBv         = @. BB*exp(-RTB/RTv)
        mFdotv      = @. -kIv*rhov*YFv - kBv*rhov^2*YFv*YXv
        mXdotv      = @. kIv*rhov*YFv + kBv*rhov^2*YFv*YXv - BR*rhov*YXv
        mPdotv      = @. BR*rhov*YXv
        
        fig = figure()
        suptitle("\$M_0 = 1.0 \\times M_{CJ} = $(round(M0, digits=3))\$")
        ax1 = subplot(1, 3, 1)
        plot(-xv, rhov ./rho1, "b", linewidth=2.5)
        plot(-xv, uv ./u1, "y", linewidth=2.5)
        plot(-xv, pv ./p1, "g", linewidth=2.5)
        plot(-xv, RTv ./RT1, "r", linewidth=2.5)
        axhline(Tcross/RT1, color="r", linestyle="--", linewidth=1.5)
        plot(-xv, Mv, "k", linewidth=2.5)
        legend([latexstring("\\rho/\\rho_1"), 
                latexstring("u/u_1"), 
                latexstring("p/p_1"),
                latexstring("T/T_1"), 
                latexstring("(RT)_{cross}/(RT)_1"), 
                latexstring("M") ], 
                loc="best")
        xlabel("x") 
        grid("on")
        
        ax2 = subplot(1, 3, 2)
        ax2.sharex(ax1)
        plot(-xv, YFv, color="b", linewidth=2.5)
        plot(-xv, YXv, color="orange", linewidth=2.5)
        plot(-xv, 1.0.-YFv.-YXv, color="g", linewidth=2.5)
        legend([latexstring("Y_{F}"), 
                latexstring("Y_{X}"),
                latexstring("Y_{P}")], 
                loc="best")
        xlabel("x") 
        grid("on")
        
        ax3 = subplot(1, 3, 3)
        ax3.sharex(ax1)
        plot(-xv, -mFdotv, color="b", linewidth=2.5)
        plot(-xv, mXdotv, color="orange", linewidth=2.5)
        plot(-xv, mPdotv, color="g", linewidth=2.5)
        legend([latexstring("-\\dot{m}_{F}"), 
                latexstring("\\dot{m}_{X}"),
                latexstring("\\dot{m}_{P}")], 
                loc="best")
        xlabel("x") 
        grid("on")
    
    end
    
    return xv, rhov, uv, pv, YFv, YXv
    
end


"""
    Overdriven(M0, gamma, q) -> M2, rho0, u0, p0, rho1, u1, p1, rho2, u2, p2

States of an overdriven detonation with Mach number `M0` (>= CJ Mach number), with the same
state numbering as [`ChapmanJouguet`](@ref). `M2` is the (subsonic) Mach number of the
burnt gas. Throws an error if `M0` is below the CJ Mach number.
"""
function Overdriven(M0::Float64, gamma::Float64, q::Float64)

    #We assume the state at the end p0=1, rho0=1, u0=M0*sqrt(gamma), YF2=0.
    p0      = 1.0
    rho0    = 1.0
    u0      = M0*sqrt(gamma)
    
    if M0 < M0_ChapmanJouguet(gamma, q)
        error("Initial Mach M0 = $M0 is lower than the minimum CJ Mach for q = $q.")
    end

    #Then, we compute the state "2" by imposing conservation equations:
    #   rho0*u0                                                 = M
    #   p0+rho0*u0^2                                            = P
    #   (gamma/(gamma-1)*p0 + 0.5*rho0*u0^2 + rho0*YF0*Q0)*u0   = HU
    #Note that rho0*Q0*YF0 = q*gamma/(gamma-1)*p0, so the last equation is
    #   ((1+q)*gamma/(gamma-1)*p0 + 0.5*rho0*u0^2)*u0           = HU
    #If we write the equations in terms of u0:
    #A*u2^2 + B*u2 + HU = 0
    M       = rho0*u0
    P       = p0+rho0*u0^2
    HU      = ((1+q)*gamma/(gamma-1.0)*p0 + 0.5*rho0*u0^2)*u0
    #
    A       = (0.5-gamma/(gamma-1.0))*M # A<0
    B       = gamma/(gamma-1.0)*P
    C       = -HU
    #

    # Physical feasibility check (Mach >= Mach_CJ)
    disc = B^2 - 4.0 * A * C
    if disc < 0.0
        error("Initial Mach M0 = $M0 is lower than the minimum CJ Mach for q = $q.")
    end

    u2      = (-B+sqrt(B^2-4*A*C))/(2*A) # We take the positive root to take the smaller u2 (overdriven detonation)
    rho2    = M/u2
    p2      = P-M*u2
    M2      = u2/sqrt(gamma*p2/rho2)
   

    #Rescale solution with entry conditions:
    #rho2    *= 1.0/rho0
    #p2      *= 1.0/p0
    #u2      = sqrt(gamma*p2/rho2)
    #
    #rho0    = 1.0
    #p0      = 1.0
    #u0      = M0*sqrt(gamma)
    
    #State at region 1:
    rho1            = rho0 / ( (gamma-1.0)/(gamma+1.0) + 2.0/(gamma+1.0)/M0^2 )
    u1              = rho0*u0/rho1
    p1              = p0+rho0*u0^2-rho1*u1^2
    
    #Check results:
#     display(rho0*u0 - rho2*u2)
#     display(p0+rho0*u0^2 - (p2+rho2*u2^2))
#     display(((1+q)*gamma/(gamma-1)*p0 + 0.5*rho0*u0^2)*u0 - 
#             (gamma/(gamma-1)*p2 + 0.5*rho2*u2^2)*u2)
#     display(rho1*u1 - rho2*u2)
#     display(p1+rho1*u1^2 - (p2+rho2*u2^2))
#     display(((1+q)*gamma/(gamma-1)*p1 + 0.5*rho1*u1^2)*u1 - 
#             (gamma/(gamma-1)*p2 + 0.5*rho2*u2^2)*u2)
    
    
return M2, rho0, u0, p0, rho1, u1, p1, rho2, u2, p2

end


"""
    FXP_Overdriven(; gamma, q, M0, RThatI, RThatB, BIhat, BBhat, Deltax, xf, YF0, YX0, PlotRes)
        -> xv, rhov, uv, pv, YFv, YXv

Same model as [`FXP_ChapmanJouguet`](@ref) (F -> X -> P kinetics) but for an overdriven
detonation of Mach number `M0` (>= CJ Mach number). The post-shock state is obtained from
the Rankine-Hugoniot relations for `M0`, whereas the activation temperatures and the
branching factor are normalized with the CJ post-shock state (RT1_CJ, rho1_CJ), so that
the kinetics are the same as in the CJ case.
"""
function FXP_Overdriven(; gamma::Float64=1.4, q::Float64=5.0, M0::Float64=6.0,
    RThatI::Float64=10.0, RThatB::Float64=10.0, #RThat=RT_a/RT_1, #1: means post-shock condition
    BIhat::Float64=2.5e-3, BBhat::Float64=4.4e5, 
    Deltax::Float64=1e-3, xf::Float64=30.0, YF0::Float64=1.0-1e-6, YX0::Float64=1e-6, 
    PlotRes::Bool=false)

    # Compute CJ reference values:
    M0_CJ, rho0, u0_CJ, p0, rho1_CJ, u1_CJ, p1_CJ, _, _, _ = ChapmanJouguet(gamma, q)

    if M0_CJ > M0
        error("Initial Mach M0 = $M0 is lower than the minimum CJ Mach for q = $q.")
    end

    u0 = M0 * sqrt(gamma * p0 / rho0)
    #Compute detonation states:
    rho1 = rho0 * ((gamma + 1.0) * M0^2) / (2.0 + (gamma - 1.0) * M0^2)
    u1   = u0 * (rho0 / rho1)
    p1   = p0 * (1.0 + (2.0 * gamma / (gamma + 1.0)) * (M0^2 - 1.0))
    
    #Characteristic variables:
#     rhoc            = rho1
#     uc              = u1
#     pc              = rhoc*uc^2
#     tc              = 1.0
#     lc              = uc*tc
    
    #Scale pre and post shock conditions:
#     rho0            /= rhoc
#     u0              /= uc
#     p0              /= pc
#     rho1            /= rhoc
#     u1              /= uc
#     p1              /= pc
#     rho2            /= rhoc
#     u2              /= uc
#     p2              /= pc
    
    #Compute RT:
    RT0             = p0/rho0
    RT1             = p1/rho1
    
    #Compute constants:
    Q               = gamma/(gamma-1)*RT0/YF0*q     #F and X formation mass enthalpy 
    RT1_CJ = p1_CJ/rho1_CJ
    RTI_CJ             = RThatI*RT1_CJ                    #Activation temperature of initiation (CJ-normalized)
    RTB_CJ             = RThatB*RT1_CJ                    #Activation temperature of chain-branching (CJ-normalized)
    tc              = 1.0                           #Characteristic time
    BR              = 1.0/tc                        #t_c = B_R^{-1} = 1.0
    BI              = BIhat * BR
    BB              = BBhat * BR/rho1_CJ
                
    Tcross = compute_FXP_Tcross(gamma=gamma, q=q, RThatB=RThatB, BBhat=BBhat, YF0=YF0,M0=M0)

    #Function to compute rho, u, p, YF, YX from vector 
    #   F   =   [   rhoY_F u, 
    #               rhoY_X u, 
    #               rho u, 
    #               p + rho u^2, 
    #               ( gamma/(gamma-1)*p + rho (Y_F+Y_X) Q + 1/2 rho u^2 ) u ]
    function ExtractVars(Fv::Vector{Float64})
    
        #YF = rho YF u / (rho u):
        YF          = Fv[1]/Fv[3]
        YX          = Fv[2]/Fv[3]
        YP          = 1-YF-YX
        
        #F3 = rho u
        #F4 = p + rho u^2 
        #F5 = ( gamma/(gamma-1)*p + 1/2 rho u^2 + rho Q (YF+YX) ) u
        #   = gamma/(gamma-1)*p*u + 1/2 F3 u^2 + (F1+F2) Q
        #Combining these equations:
        #(gamma/(gamma-1)-1/2) F3 U^2 - gamma/(gamma-1) F4 U + F5 - (F1+F2)*Q = 0
        A       = (gamma/(gamma-1)-0.5)*Fv[3]
        B       = -gamma/(gamma-1)*Fv[4]
        C       = Fv[5]-(Fv[1]+Fv[2])*Q
        if B^2-4*A*C<0
#             return NaN, NaN, NaN, NaN, NaN
        end
        u       = (-B-sqrt(B^2-4*A*C))/(2*A)    #Smaller root: subsonic branch (post-shock state); the larger root is the supersonic state 0
        
        #Therefore:
        rho         = Fv[3]/u
        p           = Fv[4]-Fv[3]*u
        
        return rho, u, p, YF, YX
        
    end
    
    #Function to march flux variables in space:
    function dF_dx(Fv::Vector{Float64}, ComputeJ::Bool)
    
        rho, u, p, YF, YX   = ExtractVars(Fv)
        RT                  = p/rho
        kI                  = BI*exp(-RTI_CJ/RT)
        kB                  = BB*exp(-RTB_CJ/RT)
        kR                  = BR
        mFdot               = - kI*rho*YF - kB*rho^2*YF*YX
        mXdot               =   kI*rho*YF + kB*rho^2*YF*YX - kR*rho*YX
        
        f                   = [ mFdot; mXdot; 0; 0; 0 ]
        J                   = zeros(0,0)
        if ComputeJ
            nVars           = length(Fv)
            J               = zeros(nVars, nVars)
            delta           = 1e-5
            for JJ=1:nVars
                Fpert       = copy(Fv)
                Fpert[JJ]   += delta
                fpert       = dF_dx(Fpert, false)[1]
                J[:,JJ]     = (fpert-f)/delta 
            end
        end
        
        return f, J
        
    end
    
    #Initialize variables:
    x_n             = 0.0
    F_n             = [ rho1*YF0*u1 ;
                        rho1*YX0*u1 ;
                        rho1*u1;
                        p1 + rho1*u1^2;
                        (gamma/(gamma-1.0)*p1 + rho1*(YF0+YX0)*Q + 0.5*rho1*u1^2)*u1 ]
                        
#     display((rho0, u0, p0, YF0, YX0))
#     display((rho1, u1, p1, YF0, YX0))
#     display(ExtractVars(F_n))
#     return
    
    #Apply RK method:
    RK              = RK_Coefficients("Ascher3")
    NF              = length(F_n)
    fm              = zeros(NF, RK.stages)
    JFm             = zeros(NF, RK.stages)
    F_np1           = zeros(NF)
    xv              = [ -xf*0.1, x_n, x_n ]
    rhov            = [ rho0, rho0, rho1 ]
    uv              = [ u0, u0, u1 ]
    pv              = [ p0, p0, p1 ]
    YFv             = [ YF0, YF0, YF0 ]
    YXv             = [ YX0, YX0, YX0 ]
        
    while x_n < xf
    
        #println("x_n=", x_n)
        
        #Compute next stage:
        fm[:,1], J      = dF_dx(F_n, true)
        JFm[:,1]        = J*F_n
        A_F             = factorize(eye(5)-Deltax*RK.AI[2,2]*J)
        for ii=2:RK.stages
            S           = zeros(NF)
            for jj=1:ii-1
                S       += RK.AE[ii,jj]*fm[:,jj] + (RK.AI[ii,jj]-RK.AE[ii,jj])*JFm[:,jj]
            end
            F_np1       .= A_F\(F_n + Deltax*S)
            fm[:,ii]    = dF_dx(F_np1, false)[1]
            JFm[:,ii]   = J*F_np1
        end
        
        #Update solution:
        x_n             += Deltax
        F_n             .= F_np1
        
        #Save solution:
        push!(xv, x_n)
        rho_n, u_n, p_n, YF_n, YX_n     = ExtractVars(F_n)
        push!(rhov, rho_n)
        push!(uv, u_n)
        push!(pv, p_n)
        push!(YFv, YF_n)
        push!(YXv, YX_n)
        
    end
    
    #Compute Mach number:
    Mv      = @. uv/sqrt(gamma*pv/rhov)
    
    #Plot each variable:
    if PlotRes
    
        Nx      = length(xv)
        #println("Exact solution: rho2=$rho2, u2=$u2, p2=$p2, YF2=0.0, M2=1.0")
        println("Numerical solution: rho2=$(rhov[Nx]), u2=$(uv[Nx]), p2=$(pv[Nx]), ", 
                "YF2=$(YFv[Nx]), M2=$(Mv[Nx])")
        #=
        figure()
#         plot(xv/lc, rhov/rho1, "b")
#         plot(xv/lc, uv/u1, "g")
#         plot(xv/lc, pv/p1, "c")
#         plot(xv/lc, (pv./rhov)/RT1, "r")
        plot(xv, rhov, "b", linewidth=2.5)
        plot(xv, uv, "g", linewidth=2.5)
        plot(xv, pv, "c", linewidth=2.5)
        plot(xv, pv./rhov, "r", linewidth=2.5)
        plot(xv, Mv, "k", linewidth=2.5)
        legend([latexstring("\\rho"), 
                latexstring("u"), 
                latexstring("p"),
                latexstring("T"), 
                latexstring("M") ], 
                loc="best")
        xlabel("x") 
        grid("on")
        
        figure()
        plot(xv, YFv, color="b", linewidth=2.5)
        plot(xv, YXv, color="orange", linewidth=2.5)
        plot(xv, 1.0.-YFv.-YXv, color="g", linewidth=2.5)
        legend([latexstring("Y_{F}"), 
                latexstring("Y_{X}"),
                latexstring("Y_{P}")], 
                loc="best")
        xlabel("x") 
        grid("on")

        # Adim with respect to post-shock conditions:
        figure()
        plot(-xv, rhov ./rho1, "b", linewidth=2.5)
        plot(-xv, uv ./u1, "y", linewidth=2.5)
        plot(-xv, pv ./p1, "g", linewidth=2.5)
        plot(-xv, pv./rhov ./RT1, "r", linewidth=2.5)
        plot(-xv, Mv, "k", linewidth=2.5)
        legend([latexstring("\\rho/\\rho_1"), 
                latexstring("u/u_1"), 
                latexstring("p/p_1"),
                latexstring("T/T_1"), 
                latexstring("M") ], 
                loc="best")
        xlabel("x") 
        grid("on")

        # Mass fractions vs -x
        figure()
        plot(-xv, YFv, color="b", linewidth=2.5)
        plot(-xv, YXv, color="orange", linewidth=2.5)
        plot(-xv, 1.0.-YFv.-YXv, color="g", linewidth=2.5)
        legend([latexstring("Y_{F}"), 
                latexstring("Y_{X}"),
                latexstring("Y_{P}")], 
                loc="best")
        xlabel("x") 
        grid("on")
        =#

        # Combined figure with shared x-axis: adim. variables, mass fractions and mdot
        RTv         = pv./rhov
        kIv         = @. BI*exp(-RTI_CJ/RTv)
        kBv         = @. BB*exp(-RTB_CJ/RTv)
        mFdotv      = @. -kIv*rhov*YFv - kBv*rhov^2*YFv*YXv
        mXdotv      = @. kIv*rhov*YFv + kBv*rhov^2*YFv*YXv - BR*rhov*YXv
        mPdotv      = @. BR*rhov*YXv
        
        fig = figure()
        suptitle("\$M_0 = $(round(M0/M0_CJ, digits=3)) \\times M_{CJ} = $(round(M0, digits=3))\$")
        ax1 = subplot(1, 3, 1)
        plot(-xv, rhov ./rho1, "b", linewidth=2.5)
        plot(-xv, uv ./u1, "y", linewidth=2.5)
        plot(-xv, pv ./p1, "g", linewidth=2.5)
        plot(-xv, RTv ./RT1, "r", linewidth=2.5)
        axhline(Tcross/RT1, color="r", linestyle="--", linewidth=1.5)
        plot(-xv, Mv, "k", linewidth=2.5)
        legend([latexstring("\\rho/\\rho_1"), 
                latexstring("u/u_1"), 
                latexstring("p/p_1"),
                latexstring("T/T_1"), 
                latexstring("(RT)_{cross}/(RT)_1"), 
                latexstring("M") ], 
                loc="best")
        xlabel("x")
        grid("on")
        #setp(ax1.get_xticklabels(), visible=false)
        
        ax2 = subplot(1, 3, 2)
        ax2.sharex(ax1)
        plot(-xv, YFv, color="b", linewidth=2.5)
        plot(-xv, YXv, color="orange", linewidth=2.5)
        plot(-xv, 1.0.-YFv.-YXv, color="g", linewidth=2.5)
        legend([latexstring("Y_{F}"), 
                latexstring("Y_{X}"),
                latexstring("Y_{P}")], 
                loc="best")
        xlabel("x")
        grid("on")
        #setp(ax2.get_xticklabels(), visible=false)
        
        ax3 = subplot(1, 3, 3)
        ax3.sharex(ax1)
        plot(-xv, -mFdotv, color="b", linewidth=2.5)
        plot(-xv, mXdotv, color="orange", linewidth=2.5)
        plot(-xv, mPdotv, color="g", linewidth=2.5)
        legend([latexstring("-\\dot{m}_{F}"), 
                latexstring("\\dot{m}_{X}"),
                latexstring("\\dot{m}_{P}")], 
                loc="best")
        #setp(ax3.get_xticklabels(), visible=false)
        xlabel("x") 
        grid("on")
    
    end
    
    return xv, rhov, uv, pv, YFv, YXv
    
end



"""
    FXP_Overdriven_Adaptive(; gamma, q, M0, RThatI, RThatB, BIhat, BBhat,
        Deltax, MaxStep, RelTol, AbsTol, xf, YF0, YX0, PlotRes)
        -> xv, rhov, uv, pv, YFv, YXv

Adaptive-step counterpart of [`FXP_Overdriven`](@ref). `Deltax` is the requested
output spacing; internal integration steps are selected using `RelTol` and `AbsTol`,
and are capped by `MaxStep`. Thus a fine output grid does not force the integrator
to take a step at every output point.
"""
function FXP_Overdriven_Adaptive(; gamma::Float64=1.4, q::Float64=5.0, M0::Float64=6.0,
    RThatI::Float64=10.0, RThatB::Float64=10.0,
    BIhat::Float64=2.5e-3, BBhat::Float64=4.4e5,
    Deltax::Float64=1e-3, MaxStep::Float64=1e-2,
    RelTol::Float64=1e-6, AbsTol::Float64=1e-9,
    xf::Float64=30.0, YF0::Float64=1.0-1e-6, YX0::Float64=1e-6,
    PlotRes::Bool=false)

    Deltax > 0.0 || throw(ArgumentError("Deltax must be positive."))
    MaxStep > 0.0 || throw(ArgumentError("MaxStep must be positive."))
    RelTol > 0.0 || throw(ArgumentError("RelTol must be positive."))
    AbsTol > 0.0 || throw(ArgumentError("AbsTol must be positive."))
    xf > 0.0 || throw(ArgumentError("xf must be positive."))

    M0_CJ, rho0, _, p0, rho1_CJ, _, p1_CJ, _, _, _ = ChapmanJouguet(gamma, q)
    M0_CJ <= M0 ||
        error("Initial Mach M0 = $M0 is lower than the minimum CJ Mach for q = $q.")

    u0 = M0 * sqrt(gamma * p0 / rho0)
    rho1 = rho0 * ((gamma + 1.0) * M0^2) / (2.0 + (gamma - 1.0) * M0^2)
    u1 = u0 * (rho0 / rho1)
    p1 = p0 * (1.0 + (2.0 * gamma / (gamma + 1.0)) * (M0^2 - 1.0))

    RT0 = p0 / rho0
    RT1 = p1 / rho1
    Q = gamma / (gamma - 1.0) * RT0 / YF0 * q
    RT1_CJ = p1_CJ / rho1_CJ
    RTI_CJ = RThatI * RT1_CJ
    RTB_CJ = RThatB * RT1_CJ
    BR = 1.0
    BI = BIhat * BR
    BB = BBhat * BR / rho1_CJ
    Tcross = compute_FXP_Tcross(
        gamma=gamma, q=q, RThatB=RThatB, BBhat=BBhat, YF0=YF0, M0=M0)

    function ExtractVars(Fv::Vector{Float64})
        A = (gamma / (gamma - 1.0) - 0.5) * Fv[3]
        B = -gamma / (gamma - 1.0) * Fv[4]
        C = Fv[5] - (Fv[1] + Fv[2]) * Q
        discriminant = B^2 - 4.0 * A * C
        discriminant >= 0.0 ||
            throw(DomainError(discriminant, "Nonphysical state in adaptive FXP integration."))
        u = (-B - sqrt(discriminant)) / (2.0 * A)
        rho = Fv[3] / u
        p = Fv[4] - Fv[3] * u
        YF = Fv[1] / Fv[3]
        YX = Fv[2] / Fv[3]
        all(isfinite, (rho, u, p, YF, YX)) && rho > 0.0 && p > 0.0 ||
            throw(DomainError(Fv, "Nonphysical state in adaptive FXP integration."))
        return rho, u, p, YF, YX
    end

    function dF_dx(Fv::Vector{Float64}, ComputeJ::Bool)
        rho, _, p, YF, YX = ExtractVars(Fv)
        RT = p / rho
        kI = BI * exp(-RTI_CJ / RT)
        kB = BB * exp(-RTB_CJ / RT)
        mFdot = -kI * rho * YF - kB * rho^2 * YF * YX
        mXdot = kI * rho * YF + kB * rho^2 * YF * YX - BR * rho * YX
        f = [mFdot, mXdot, 0.0, 0.0, 0.0]
        ComputeJ || return f, zeros(0, 0)

        J = zeros(5, 5)
        delta = 1e-5
        for JJ in eachindex(Fv)
            Fpert = copy(Fv)
            Fpert[JJ] += delta
            fpert, _ = dF_dx(Fpert, false)
            J[:, JJ] = (fpert - f) / delta
        end
        return f, J
    end

    RK = RK_Coefficients("Ascher3")
    function rk_step(F::Vector{Float64}, h::Float64)
        fm = zeros(5, RK.stages)
        JFm = zeros(5, RK.stages)
        Fnext = zeros(5)
        fm[:, 1], J = dF_dx(F, true)
        JFm[:, 1] = J * F
        A_F = factorize(eye(5) - h * RK.AI[2, 2] * J)
        for ii in 2:RK.stages
            S = zeros(5)
            for jj in 1:ii-1
                S .+= RK.AE[ii, jj] .* fm[:, jj] .+
                      (RK.AI[ii, jj] - RK.AE[ii, jj]) .* JFm[:, jj]
            end
            Fnext .= A_F \ (F + h * S)
            fm[:, ii], _ = dF_dx(Fnext, false)
            JFm[:, ii] = J * Fnext
        end
        return Fnext
    end

    function primitive_state(F::Vector{Float64})
        rho, u, p, YF, YX = ExtractVars(F)
        return [rho, u, p, YF, YX]
    end

    x = 0.0
    F = [rho1 * YF0 * u1,
         rho1 * YX0 * u1,
         rho1 * u1,
         p1 + rho1 * u1^2,
         (gamma / (gamma - 1.0) * p1 + rho1 * (YF0 + YX0) * Q +
          0.5 * rho1 * u1^2) * u1]
    xv = [-0.1 * xf, 0.0, 0.0]
    rhov = [rho0, rho0, rho1]
    uv = [u0, u0, u1]
    pv = [p0, p0, p1]
    YFv = [YF0, YF0, YF0]
    YXv = [YX0, YX0, YX0]

    output_index = 1
    output_x = min(Deltax, xf)
    h = min(MaxStep, max(xf / 100.0, Deltax))
    minimum_step = 16.0 * eps(Float64) * max(xf, 1.0)
    while x < xf
        h = min(h, MaxStep, xf - x)
        h >= minimum_step ||
            error("Adaptive FXP integration step underflow at x=$x; relax tolerances or reduce MaxStep.")

        accepted = false
        Fnext = F
        error_norm = Inf
        try
            Ffull = rk_step(F, h)
            Fhalf = rk_step(F, h / 2.0)
            Fnext = rk_step(Fhalf, h / 2.0)

            state_full = primitive_state(Ffull)
            state_next = primitive_state(Fnext)
            state_old = primitive_state(F)
            error_norm = maximum(
                abs(state_next[ii] - state_full[ii]) / 7.0 /
                (AbsTol + RelTol * max(abs(state_old[ii]), abs(state_next[ii])))
                for ii in eachindex(state_next))
            accepted = isfinite(error_norm) && error_norm <= 1.0
        catch err
            if !(err isa DomainError)
                rethrow()
            end
        end

        if accepted
            xnext = x + h
            fstart, _ = dF_dx(F, false)
            fend, _ = dF_dx(Fnext, false)
            while output_index * Deltax <= xnext + 8.0 * eps(max(xnext, 1.0)) ||
                  (xnext == xf && output_x <= xf)
                output_x = min(output_index * Deltax, xf)
                theta = clamp((output_x - x) / h, 0.0, 1.0)
                theta2 = theta^2
                theta3 = theta^3
                Fout = (2.0 * theta3 - 3.0 * theta2 + 1.0) .* F .+
                       (theta3 - 2.0 * theta2 + theta) * h .* fstart .+
                       (-2.0 * theta3 + 3.0 * theta2) .* Fnext .+
                       (theta3 - theta2) * h .* fend
                rho_out, u_out, p_out, YF_out, YX_out = ExtractVars(Fout)
                push!(xv, output_x)
                push!(rhov, rho_out)
                push!(uv, u_out)
                push!(pv, p_out)
                push!(YFv, YF_out)
                push!(YXv, YX_out)
                output_index += 1
                if output_x >= xf
                    break
                end
                output_x = min(output_index * Deltax, xf)
            end
            F = Fnext
            x = xnext
        end

        factor = accepted && error_norm > 0.0 ?
            clamp(0.9 * error_norm^(-0.25), 0.2, 2.0) :
            (accepted ? 2.0 : 0.25)
        h *= factor
    end

    if PlotRes
        Mv = @. uv / sqrt(gamma * pv / rhov)
        RTv = pv ./ rhov
        kIv = @. BI * exp(-RTI_CJ / RTv)
        kBv = @. BB * exp(-RTB_CJ / RTv)
        mFdotv = @. -kIv * rhov * YFv - kBv * rhov^2 * YFv * YXv
        mXdotv = @. kIv * rhov * YFv + kBv * rhov^2 * YFv * YXv - BR * rhov * YXv
        mPdotv = @. BR * rhov * YXv

        figure()
        suptitle("\$M_0 = $(round(M0 / M0_CJ, digits=3)) \\times M_{CJ} = $(round(M0, digits=3))\$")
        ax1 = subplot(1, 3, 1)
        plot(-xv, rhov ./ rho1, "b")
        plot(-xv, uv ./ u1, "y")
        plot(-xv, pv ./ p1, "g")
        plot(-xv, RTv ./ RT1, "r")
        axhline(Tcross / RT1, color="r", linestyle="--")
        plot(-xv, Mv, "k")
        xlabel("x")
        grid("on")
        ax2 = subplot(1, 3, 2)
        ax2.sharex(ax1)
        plot(-xv, YFv, color="b")
        plot(-xv, YXv, color="orange")
        plot(-xv, 1.0 .- YFv .- YXv, color="g")
        xlabel("x")
        grid("on")
        ax3 = subplot(1, 3, 3)
        ax3.sharex(ax1)
        plot(-xv, -mFdotv, color="b")
        plot(-xv, mXdotv, color="orange")
        plot(-xv, mPdotv, color="g")
        xlabel("x")
        grid("on")
    end

    return xv, rhov, uv, pv, YFv, YXv
end


"""
    compute_FXP_Tcross(; gamma=1.4, q=5.0, RThatB=10.0, BBhat=4.4e5, YF0=1.0, M0=nothing)

Computes the crossover temperature `RT_cross` (below which the chain-branching rate
is lower than the recombination rate) for both CJ and overdriven (OD) detonations.
If `M0` is `nothing` or equals `M0_CJ`, the CJ state is assumed.
"""
function compute_FXP_Tcross(; gamma::Float64=1.4, q::Float64=5.0, RThatB::Float64=10.0, 
    BBhat::Float64=4.4e5, YF0::Float64=1.0, M0::Union{Float64, Nothing}=nothing)

    # Reference CJ state
    M0_CJ, rho0, u0_CJ, p0, rho1_CJ, _, p1_CJ, _, _, _ = ChapmanJouguet(gamma, q)
    RT1_CJ = p1_CJ / rho1_CJ
    RTB_CJ = RThatB * RT1_CJ

    # Evaluate post-shock density and temperature for the requested Mach
    if isnothing(M0) || M0 ≈ M0_CJ
        rho1 = rho1_CJ
        RT1  = RT1_CJ
    else
        # Permitir cualquier choque supersónico (incluso choque precursor pre-DDT)
        if M0 <= 1.0
            error("Mach number M0 = $M0 must be supersonic (M0 > 1.0)")
        end
        u0   = M0 * sqrt(gamma * p0 / rho0)
        rho1 = rho0 * ((gamma + 1.0) * M0^2) / (2.0 + (gamma - 1.0) * M0^2)
        p1   = p0 * (1.0 + (2.0 * gamma / (gamma + 1.0)) * (M0^2 - 1.0))
        RT1  = p1 / rho1
    end

    # Branching / Recombination equilibrium: k_B * rho1 * YF0 = k_R
    arg = BBhat * YF0 * (rho1 / rho1_CJ)
    if arg <= 1.0
        error("Pre-exponential factor BBhat is too low to sustain chain-branching at this shock strength.")
    end

    RT_cross = RTB_CJ / log(arg)

    return RT_cross
end