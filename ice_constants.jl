#=
Structures to contain constants and other useful properties for ice sheet modeling
=# 

using Parameters

@with_kw mutable struct Material
    n::Int = 3
    B::Float64 = 3.2e8
    ρ_w::Float64 = 1030.0
    ρ_i::Float64 = 916.0
    g::Float64 = 9.81 
    seconds_in_year::Float64 = 365.2422*24*60*60
    di::Float64 = (1-ρ_i/ρ_w)
end

@with_kw mutable struct Shelf 
    velocity 
    stress 
    thick
    thick_cellCenter 
end

mutable struct Stress
    xx
    yy
    xy
end

mutable struct Velocity
    u
    v
end

#mutable struct ice_thick 
#    thick
#    thick_cellCenter
#end



@with_kw struct Dimensionless
    mach_number::Float64 = 0.1
    froude_number::Float64 = 1.0
    argand_number::Float64 = 1.0
    mach_number_squared::Float64 = mach_number^2
    froude_number_squared::Float64 = froude_number^2
    argand_upon_mach2::Float64 = argand_number/mach_number^2
    froude_upon_mach2::Float64 = froude_number^2/mach_number^2
end


#=
function visc(τ,thick_edgeFace;n=3.0)
    local Σ = sqrt.(τ.xx.^2 .+ τ.yy.^2 .+ τ.xy.^2 .+ τ.xx.*τ.xy)
    local η = ((abs.(Σ)./thick_edgeFace)).^(1-n)
    η[:,1]=η[:,2]
    η[:,end]=η[:,end-1]
    return η
end

function visc_plas(τ,thick_edgeFace;n=3.0,eps=0.0,tau_y=200.0,tau_min=10.0,eps_crit=0.1)
    local tau_max = max.(tau_y*(1-eps./eps_crit),tau_min)
    local Σ = sqrt.(τ.xx.^2 .+ τ.yy.^2 .+ τ.xy.^2 .+ τ.xx.*τ.xy)./thick_edgeFace
    local η_visc = ((abs.(Σ)) .+ 1e-12).^(1-n)
    local η_plas = tau_max./(Σ.^n)
    η = min.(η_plas,η_visc)
    η[:,1]=η[:,2]
    η[:,end]=η[:,end-1]
    return η
end
=#








