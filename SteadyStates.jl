#include("ice_constants.jl")
function steady_state_old(x,L0,accum_rate,B0,H_input,U_input)
    A0 = (1/B0)^n;                      # Rate constant for ice
    a0 = accum_rate 
    C = (rho_ice*g*(1-rho_ice/rho_w)/4/B0)^n;
    if a0 == 0
        h = ((n+1)*C./(H_input)./(U_input)*L0*(x) .+ (H_input)^(-n-1)).^(-1/(n+1));
        u = (H_input*U_input)./h;
    elseif a0>0
        M =  a0;
        h = (C./M .- (U_input).^(n+1).*(C./M.*(H_input).^(n+1)-1)./(L0*x*M .+H_input*U_input).^(n+1)).^(-1/(n+1));
        u = (a0*L0*x .+ H_input*U_input)./h;
    elseif a0<0
        M = -a0;
        h = (-C/M .+ (U_input).^(n+1)*(1 .+ C./M*(H_input).^(n+1))./(H_input*U0 .- L0*x*M).^(n+1)).^(-1/(n+1));
        u = (a0*L0*x .+ H_input*U_input)./h;
    end
    Exx = C*h.^n;
    S = B0.*Exx.^(1/3)./(rho_ice*g*h)*4/(1-rho_ice/rho_w);
    
    return h,u,Exx,S
end

function steady_state(x,L0,accum_rate,H_input,U_input,mat)
    #=
    Function to calculate analytic ice shelf profile, velocity and strain rate
    Inputs: x - horizontal position vector
            L0 - length of ice shelf
            accum rate - accumulation rate (m/a)
            B0 - rate constant (s^{1/3} Pa)
            H_input - ice thickness at the grounding line (m)
            U_input - ice velocity at the grounding line (m/a)
    Returns:
            h - ice thickness evaluated at x (m)
            u - velocity evaluated at x (m/a)
            Exx - strain rate evaluated at x (1/a)
            S - deviatoric stress evaluated at x (Pa)
    =#


    n = mat.n
    B0 = mat.B
    rho_w = mat.ρ_w
    rho_ice = mat.ρ_i
    g = mat.g
    A0 = (1/B0)^n;                      # Rate constant for ice
    a0 = accum_rate 
    C = (rho_ice*g*(1-rho_ice/rho_w)/4/B0)^n;
    if a0 == 0
        h = ((n+1)*C./(H_input)./(U_input)*L0*(x) .+ (H_input)^(-n-1)).^(-1/(n+1));
        u = (H_input*U_input)./h;
    elseif a0>0
        M =  a0;
        h = (C./M .- (U_input).^(n+1).*(C./M.*(H_input).^(n+1)-1)./(L0*x*M .+H_input*U_input).^(n+1)).^(-1/(n+1));
        u = (a0*L0*x .+ H_input*U_input)./h;
    elseif a0<0
        M = -a0;
        h = (-C/M .+ (U_input).^(n+1)*(1 .+ C./M*(H_input).^(n+1))./(H_input*U_input .- L0*x*M).^(n+1)).^(-1/(n+1));
        u = (a0*L0*x .+ H_input*U_input)./h;
    end
    Exx = C*h.^n;
    S = B0.*Exx.^(1/3)./(rho_ice*g*h)*4/(1-rho_ice/rho_w);
    
    return h,u,Exx,S
end

function steady_state(x,y,L0,accum_rate,H_input,U_input,mat)
    #=
    Function to calculate analytic ice shelf profile, velocity and strain rate
    Inputs: x - horizontal position vector
            y - horizontal position vector
            L0 - length of ice shelf
            accum rate - accumulation rate (m/a)
            B0 - rate constant (s^{1/3} Pa)
            H_input - ice thickness at the grounding line (m)
            U_input - ice velocity at the grounding line (m/a)
    Returns:
            h - ice thickness evaluated at x and y (m)
            u - velocity evaluated at x at x and y(m/a)
    =#

    ha,ua,Exx=steady_state(x,L0,accum_rate,H_input,U_input,mat)

    # Blow up the 1D solution to 2D
    thick = [ha[i] for j in 1:length(y), i in 1:length(x)]
    u = [ua[i] for j in 1:length(y), i in 1:length(x)]
    v = zeros(size(u))
    return thick,u,v
end


