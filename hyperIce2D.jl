#=
Modules to set boundary conditions and update velocity and stresses using Lax Wendroff
=#
using Parameters
import Statistics

@with_kw mutable struct visc_prms 
    n::Int = 3 # Flow law exponent
    eps_crit::Float64 = 0.1 # Critical strain
    tau_y::Float64 = 1.0 # Non-dimensional yield strength
    tau_min::Float64 = 0.1 # Non-dimensional minimum yield strength
    μ::Float64 = 0.0  # Friction coefficent
    alpha::Float64 = 0.0 # Relative importance of largest principal stress vs effective stress
    eps
end

function visc(τ,thick_edgeFace;n=3.0)
    local Σ = sqrt.(τ.xx.^2 .+ τ.yy.^2 .+ τ.xy.^2 .+ τ.xx.*τ.yy .+ 1e-12.^2)
    local η = ((abs.(Σ)./thick_edgeFace)).^(1-n)
    η[:,1]=η[:,2]
    η[:,end]=η[:,end-1]
    return η
end

function get_max_stress(τ::Stress,thick_edgeFace::AbstractArray{<:Number},visc_prms::visc_prms)
    local t11= (0.5*(τ.xx+τ.yy) + 0.5*sqrt.(τ.xx.^2 +τ.yy.^2 - 2*τ.xx.*τ.yy +4*τ.xy.^2))./thick_edgeFace
    local t22= (0.5*(τ.xx+τ.yy) - 0.5*sqrt.(τ.xx.^2 +τ.yy.^2 - 2*τ.xx.*τ.yy +4*τ.xy.^2))./thick_edgeFace
    #db = min.((2*t11+t22)*0.22233009708737866./thick_edgeFace,0.99)
    #S0 =(2*t11+t22)./4*thick_edgeFace
    #db = (2*t11+t22)./4.0./thick_edgeFace
    db = 0.0
    eps_crit = visc_prms.eps_crit.*thick_edgeFace
    tau_y = visc_prms.tau_y.*(1.0 .- db) .+ visc_prms.μ*thick_edgeFace#.*(1.0.-exp.(-thick_edgeFace))
    return tau_max = max.(tau_y.*(1.0.-visc_prms.eps./eps_crit),visc_prms.tau_min)
end

function eps_plastic(τ::Stress,thick_edgeFace::AbstractArray{<:Number},visc_prms::visc_prms)
    local Σ = sqrt.(τ.xx.^2 .+ τ.yy.^2 .+ τ.xy.^2 .+ τ.xx.*τ.yy)./thick_edgeFace
    local t11= (0.5*(τ.xx+τ.yy) + 0.5*sqrt.(τ.xx.^2 +τ.yy.^2 - 2*τ.xx.*τ.yy +4*τ.xy.^2))./thick_edgeFace
    local t22= (0.5*(τ.xx+τ.yy) - 0.5*sqrt.(τ.xx.^2 +τ.yy.^2 - 2*τ.xx.*τ.yy +4*τ.xy.^2))./thick_edgeFace
    #local T = sqrt.(0.5*t11.^2.0.*(t11.>0)+0.5*t22.^2.0.*(t22.>0))
    #local Σ = sqrt.(0.5*τ.xx.^2 .+ 0.5*τ.yy.^2 .+ τ.xy.^2 )./thick_edgeFace
    eps = Σ.^(visc_prms.n-1).*(visc_prms.alpha*(t11) .+ (1-visc_prms.alpha).*Σ)
    #eps = T.^(visc_prms.n-1).*(visc_prms.alpha*(t11) .+ (1-visc_prms.alpha).*T)
    #eps = Σ.^(visc_prms.n-1).*(t11 + 0.5*t22)
    #eps =(t11 + 0.5*t22)
    return  eps.*(eps.>0)
end


function visc(τ::Stress,thick_edgeFace::AbstractArray{<:Number},visc_prms::visc_prms)
    #eps_crit = visc_prms.eps_crit.*thick_edgeFace
    #tau_y = visc_prms.tau_y#.*(1.0.-exp.(-thick_edgeFace))
    #local tau_max = max.((tau_y).*(1.0.-visc_prms.eps./eps_crit),visc_prms.tau_min)
    local tau_max = get_max_stress(τ,thick_edgeFace,visc_prms) 
    local Σ = sqrt.(τ.xx.^2 .+ τ.yy.^2 .+ τ.xy.^2 .+ τ.xx.*τ.yy .+ 1e-12.^2)./thick_edgeFace
    #local t11= (0.5*(τ.xx+τ.yy) + 0.5*sqrt.(τ.xx.^2 +τ.yy.^2 - 2*τ.xx.*τ.yy+4*τ.xy.^2))./thick_edgeFace
    #local epsI = max.(Σ.^(visc_prms.n-1).*t11,1e-12)*0.9
    #local epsI = Σ.^(visc_prms.n).*(t11.>0.0)
    epsI = eps_plastic(τ,thick_edgeFace,visc_prms)
    local η = ((abs.(Σ))).^(1-visc_prms.n)
    local η_plas = tau_max./epsI
    η = min.(η_plas,η)
    #η[:,1]=η[:,2]
    #η[:,end]=η[:,end-1]
    return η
end

function update_strain(τ::Stress,thick_edgeFace::AbstractArray{<:Number},visc_prms::visc_prms,dt)
    local tau_max = get_max_stress(τ,thick_edgeFace,visc_prms)
    local Σ = sqrt.(τ.xx.^2 .+ τ.yy.^2 .+ τ.xy.^2 .+ τ.xx.*τ.yy)./thick_edgeFace
    local η_visc = ((abs.(Σ)) ).^(1-visc_prms.n)
    local epsII = Σ.^visc_prms.n
    #local t11= (0.5*(τ.xx+τ.yy) + 0.5*sqrt.(τ.xx.^2 +τ.yy.^2 - 2*τ.xx.*τ.yy+4*τ.xy.^2))./thick_edgeFace
    #local t22= (0.5*(τ.xx+τ.yy) - 0.5*sqrt.(τ.xx.^2 +τ.yy.^2 - 2*τ.xx.*τ.yy+4*τ.xy.^2))./thick_edgeFace
    #local eps1 = Σ.^(visc_prms.n-1.0).*t11
    #local eps2 = Σ.^(visc_prms.n-1.0).*t22
    #local eps = max.(eps1+eps2,0.0)
    epsI = eps_plastic(τ,thick_edgeFace,visc_prms)
    local η_plas = tau_max./epsI
    visc_prms.eps = min.(visc_prms.eps + epsI.*(η_plas .< η_visc)*dt,visc_prms.eps_crit)
    #local t11= (0.5*(τ.xx+τ.yy) + 0.5*sqrt.(τ.xx.^2 +τ.yy.^2 - 2*τ.xx.*τ.yy +4*τ.xy.^2))./thick_edgeFace
    #local t22= (0.5*(τ.xx+τ.yy) - 0.5*sqrt.(τ.xx.^2 +τ.yy.^2 - 2*τ.xx.*τ.yy +4*τ.xy.^2))./thick_edgeFace
    #db = min.((2*t11+t22)*0.22233009708737866./thick_edgeFace,0.99)
    #S0 = (t11+0.5*t22)./2.0./thick_edgeFace;
    #db = (2*t11+t22)./4.0./thick_edgeFace
    #visc_prms.eps[S0.>visc_prms.μ] .= visc_prms.eps_crit/100
    
    return visc_prms.eps
end

function set_bc(u,v,txx,tyy,txy,thick,U_inflow,V_inflow)


    # Extrapolate v  velocity boundary conditions
    u[1,:] = u[2,:]
    u[end,:]=u[end-1,:]
    u[:,end] = u[:,end-1]
    u[:,1]=2*U_inflow .- u[:,2] # Average across cells = grounding line v
    #u[:,1]=u[:,2]
    # Extrapolate v  velocity boundary conditions
    v[1,:] = 2*V_inflow.-v[2,:]
    v[end,:]=2*V_inflow.-v[end-1,:]
    v[:,1]= v[:,2] # Average across cells = margin vel
    v[:,end]=v[:,end-1] # Average across cells = margin vel

    # Extrapolate to apply stress boundary conditions
    txx[1,:] = txx[2,:]
    txx[end,:]=txx[end-1,:]
    h0 = 0.5*(thick[:,end]+thick[:,end-1])
    txx[:,end]= 2*h0.^2-txx[:,end-1]
    #txx[:,end]= txx[:,end-1]
    txx[:,1]= txx[:,2]

    tyy[:,1] = tyy[:,2];tyy[:,end] = tyy[:,end-1]    
    tyy[1,:] = tyy[2,:]; tyy[end,:] = tyy[end-1,:] 
    #tyy[2:end,end]= -tyy[2:end,end-1]

    txy[:,end] = txy[:,end ];txy[:,1] = txy[:,2]
    txy[1,:] = txy[2,:];txy[end,:] = txy[end-1,:]
end

function set_thick_bc(thick)
    thick[1,:] = thick[2,:];
    thick[end,:]=thick[end-1,:]
    thick[:,end] = thick[:,end-1]
    thick[:,1] = thick[:,2]
    thick[:,1]=2 .- thick[:,2]
end

function initialize2D(grid,ua,ha,H_inflow)
    xx=grid.xx 
    yy=grid.yy
    Nx=grid.Nx 
    Ny=grid.Ny

    # Expand the analytic function to 2D 
    thick_cellCenter = [ha[i] for j in 1:length(yy), i in 1:length(xx)]
    u_cellCenter = [ua[i] for j in 1:length(yy), i in 1:length(xx)]

    # Initialize variables
    v_cellCenter = zeros(size(u_cellCenter))
    Txx = thick_cellCenter.^2
    Tyy = zeros(size(Txx))
    Txy = zeros(size(Txx))
    u = zeros(Ny+2,Nx+2);v = zeros(Ny+2,Nx+2);
    txx = zeros(Ny+2,Nx+2); tyy=zeros(Ny+2,Nx+2); txy=zeros(Ny+2,Nx+2);
    thick = zeros(Ny+2,Nx+2);
    #u1 = 0.5*(u[1:end-1,:]+u[2:end,:])
    #v1 = 0.5*(v[1:end-1,:]+v[2:end,:])
    #u[2:end-1,2:end-1] = 0.5*(u1[:,1:end-1]+u1[:,2:end])
    #v[2:end-1,2:end-1] = 0.5*(v1[:,1:end-1]+v1[:,2:end])
    #thick1 = 0.5*(thick_cellCenter[1:end-1,:]+thick_cellCenter[2:end,:])
    #thick[2:end-1,2:end-1] = 0.5*(thick1[:,1:end-1]+thick1[:,2:end])
    #txx1 = 0.5*(Txx[1:end-1,:]+Txx[2:end,:])
    #txx[2:end-1,2:end-1] = 0.5*(txx1[:,1:end-1]+txx1[:,2:end])
    #tyy1 = 0.5*(Tyy[1:end-1,:]+Tyy[2:end,:])
    #tyy[2:end-1,2:end-1] = 0.5*(tyy1[:,1:end-1]+tyy1[:,2:end])
    #txy1 = 0.5*(Txy[1:end-1,:]+Txy[2:end,:])
    #txy[2:end-1,2:end-1] = 0.5*(txy1[:,1:end-1]+txy1[:,2:end])
    for i in 2:Ny+1
        for j in 2:Nx+1
            u[i,j]=0.25*(u_cellCenter[i-1,j-1]+u_cellCenter[i,j-1]+u_cellCenter[i-1,j]+u_cellCenter[i,j])
            v[i,j]=0.25*(v_cellCenter[i-1,j-1]+v_cellCenter[i,j-1]+v_cellCenter[i-1,j]+v_cellCenter[i,j])
            thick[i,j]=0.25*(thick_cellCenter[i-1,j-1]+thick_cellCenter[i,j-1]+thick_cellCenter[i-1,j]+thick_cellCenter[i,j])
            txx[i,j]=0.25*(Txx[i-1,j-1]+Txx[i,j-1]+Txx[i-1,j]+Txx[i,j])
            tyy[i,j]=0.25*(Tyy[i-1,j-1]+Tyy[i,j-1]+Tyy[i-1,j]+Tyy[i,j])
            txy[i,j]=0.25*(Txy[i-1,j-1]+Txy[i,j-1]+Txy[i-1,j]+Txy[i,j])
        end
    end
    vel = Velocity(u,v)
    #vel_cellCenter = Velocity(u_cellCenter,v_cellVenter)
    T = Stress(txx,tyy,txy)

    thick[1,:] = thick[2,:];
    thick[end,:]=thick[end-1,:]
    thick[:,end] = thick[:,end-1]
    #
    
    
    
    thick[:,1] = thick[:,2]
    thick[:,1]=2*H_inflow .- thick[:,2]

    #thick[:,end]= thick[:,end-1] .+ (thick[:,end-1]-thick[:,end-2])

    #h = ice_thick(thick,thick_cellCenter)

    return vel, T, thick,thick_cellCenter
end

function time_step(num_time_steps,T,η,thick,thick_edgeFace,U_inflow,numeric,prms)
    for timeCount in 1:num_time_steps
        # Calculate viscosity
    end
end

function cellCenter(var)
    # Calculate variable on staggered grid (sg)
    var_cellCenter =  0.5*(var[:,2:end]+var[:,1:end-1])
    var_cellCenter =  0.5*(var_cellCenter[2:end,:]+var_cellCenter[1:end-1,:])
    return var_cellCenter
end


function grad_thick(thick,grid)

   dx = grid.dx
   Nx = grid.Nx 
   Ny = grid.Ny
   # Calculate thickness gradient on staggered grid
   thick_grad_x1 = hcat(zeros(Ny,1),(thick[1:end-1,3:end] - thick[1:end-1,1:end-2]),(thick[1:end-1,end-2]-4*thick[1:end-1,end-1]+3*thick[1:end-1,end]))
   thick_grad_x2 = hcat(zeros(Ny,1),(thick[2:end,3:end] - thick[2:end,1:end-2]),(thick[2:end,end-2]-4*thick[2:end,end-1]+3*thick[2:end,end]))
   thick_grad_x = 0.5*(thick_grad_x1+thick_grad_x2)/2/dx#


   thick_grad_y1 = vcat(zeros(1,Nx),(thick[3:end,1:end-1] - thick[1:end-2,1:end-1]),zeros(1,Nx))
   thick_grad_y2 = vcat(zeros(1,Nx),(thick[3:end,2:end] - thick[1:end-2,2:end]),zeros(1,Nx))
   thick_grad_y = 0.5*(thick_grad_y1 .+ thick_grad_y2)/2/dx

   return thick_grad_x, thick_grad_y
end

function update_vel(vel,T,h,U_inflow,V_inflow,grid,prms,visc_prms;num_steps=1,mu=1.0,tau_margin=0.0,mask = false,beta=false)
    Nx = grid.Nx
    Ny = grid.Ny
    if mask == false
        mask = ones(Nx+2)
    end
    u = vel.u 
    v = vel.v 

    txx = T.xx 
    tyy = T.yy 
    txy = T.xy

    
    dt = grid.dt 
    dx = grid.dx
    froude_number_squared = prms.froude_number_squared
    argand_upon_mach2 = prms.argand_upon_mach2
    mach_number_squared = prms.mach_number_squared
    froude_upon_mach2 = prms.froude_upon_mach2
    argand_number = prms.argand_number


    # Square and ratio of some non-dimensional parameters
    #froude_number_squared = froude_number^2
    #mach_number_squared = mach_number^2
    #froude_upon_mach2 = froude_number_squared/mach_number_squared
    #argand_upon_mach2 = argand_number / mach_number_squared


    thick_edgeFace = h
    thick = cellCenter(h)
    thick_cellCenter = thick
    #thick = ice_thick.thick_cellCenter
    #thick_grad_x,thick_grad_y = grad_thick(thick_cellCenter,grid)
    
    thick_sy =  0.5*(thick_cellCenter[:,2:end]+thick_cellCenter[:,1:end-1])
    thick_sx =  0.5*(thick_cellCenter[2:end,:]+thick_cellCenter[1:end-1,:])

    h_mid_x = 0.5*(thick_edgeFace[:,1:end-1]+thick_edgeFace[:,2:end])
    h_mid_y = 0.5*(thick_edgeFace[1:end-1,:]+thick_edgeFace[2:end,:])

    #thick_sx,thick_sy = cellCenter(ice_thick.thick_cellCenter)
    #h_end=0.5*(thick_cellCenter[1:end-1,end]+thick_cellCenter[2:end,end])

    #grad_thick = [0;(thick[3:end] - thick[1:end-2]);thick[end-2]-4*thick[end-1]+3*thick[end]]/2/dx;


    #set_bc(u,v,T.xx,T.yy,T.xy,thick,U_inflow,V_inflow)


    t=0
    σxx = T.xx + 0.5*T.yy
    σyy = T.yy + 0.5*T.xx
    σxy = T.xy

    
    

    
    #mask[end-25:end].=0.0
  
    for timeCount in 1:num_steps
       
        η = mu.*visc(T,thick_edgeFace,visc_prms)
        #η[:,1]=η[:,2]
        #η[:,end]=η[:,end-1]
        #η[1,:]=η[2,:]
        #η[end,:]=η[end-1,:]     
   
        
        
        
        # Update time
       
        #println(t)
        #σxy = 0.5*txy
        #set_bc(u,v,T.xx,T.yy,T.xy,thick,U_inflow,V_inflow)

        #σyy = T.yy + 0.5*T.xx
        #σxx = T.xx + 0.5*T.yy 
       
    
        
        # Calculate viscosity
        σxx_minus_h2 = σxx- thick_edgeFace.^2
        σyy_minus_h2 = σyy- thick_edgeFace.^2
        
        σxx_minus_h2[:,1]=σxx_minus_h2[:,2]
        #σxx_minus_h2[:,end]=σxx_minus_h2[:,end-1]
        #σxx_minus_h2[1,:]=σxx_minus_h2[2,:]
        #σxx_minus_h2[end,:]=σxx_minus_h2[end-1,:]

        # Calculate diffusion terms (hu_x)_x, (hu_y)_y
        #hu_xx,hu_yy = diffusivity(thick_edgeFace,u,dx)
        # Calculate diffusion terms (hv_x)_x, (hv_y)_y
        #hv_xx,hv_yy = diffusivity(thick_edgeFace,v,dx)

       
         # Calculate [σyy_x)/thick]_x,[σyy_y)/thick]_y
         #σxx_minus_h2_xx,σxx_minus_h2_yy = diffusivity(1.0./thick_edgeFace,σxx_minus_h2,dx)
         # Calculate [σyy_x)/thick]_x,[σyy_y)/thick]_y
         #σyy_minus_h2_xx,σyy_minus_h2_yy = diffusivity(1.0./thick_edgeFace,σyy_minus_h2,dx)



        hv_y = (v[3:end,:]-v[1:end-2,:]).*thick_edgeFace[2:end-1,:]/(2*dx)
        hu_y = (u[3:end,:]-u[1:end-2,:]).*thick_edgeFace[2:end-1,:]/(2*dx)
        hu_x = (u[:,3:end]-u[:,1:end-2]).*thick_edgeFace[:,2:end-1]/(2*dx)
        hv_x = (v[:,3:end]-v[:,1:end-2]).*thick_edgeFace[:,2:end-1]/(2*dx)

        
        σxy_y   = (σxy[3:end,:]-σxy[1:end-2,:])./thick_edgeFace[2:end-1,:]/(2*dx) 
        σxy_x   = (σxy[:,3:end]-σxy[:,1:end-2])./thick_edgeFace[:,2:end-1]/(2*dx)

       

        #σxx_minus_h2_x = (σxx_minus_h2[:,2:end]-σxx_minus_h2[:,1:end-1])./(0.5*(thick_edgeFace[:,1:end-1]+thick_edgeFace[:,2:end]))
        #right_edge = (σxx_minus_h2[:,end-2]-4*σxx_minus_h2[:,end-1]+3*σxx_minus_h2[:,end])./thick_edgeFace[:,end]
        #σxx_minus_h2_x = hcat(zeros(grid.Ny+2,1),σxx_minus_h2_x,right_edge)/(dx)


        #σxx_minus_h2_x[:,1] .= 0.0 # Assume gradient vanishes at inflow boundary . . .
        #σxx_minus_h2_xx = (σxx_minus_h2_x[:,2:end]-σxx_minus_h2_x[:,1:end-1])/(dx)

        σxx_minus_h2_x = (σxx_minus_h2[:,3:end]-σxx_minus_h2[:,1:end-2])./thick_edgeFace[:,2:end-1]/(2*dx)
        σxx_minus_h2_x[:,1] .= (σxx_minus_h2[:,3]- σxx_minus_h2[:,2])/dx./thick_edgeFace[:,2]
        σxx_minus_h2_x[:,end] .= (σxx_minus_h2[:,end]- σxx_minus_h2[:,end-1])/dx./thick_edgeFace[:,end-1]


        σyy_minus_h2_x = (σyy_minus_h2[:,3:end]-σyy_minus_h2[:,1:end-2])./thick_edgeFace[:,2:end-1]/(2*dx)
        σyy_minus_h2_x[:,1] .= (σyy_minus_h2[:,3]- σyy_minus_h2[:,2])/dx./thick_edgeFace[:,2]
        σyy_minus_h2_x[:,end] .= (σyy_minus_h2[:,end]- σyy_minus_h2[:,end-1])/dx./thick_edgeFace[:,end-1]
        #σxx_minus_h2_x[:,1] = σxx_minus_h2_x[:,2]
        
        #σxx_minus_h2_x[:,1] .= 0
        #right_edge = ((σxx_minus_h2[:,end-2]-4*σxx_minus_h2[:,end-1]+3*σxx_minus_h2[:,end])./thick_edgeFace[:,end-1])/(2*dx)
        #σxx_minus_h2_x[:,end] = right_edge
        #σxx_minus_h2_x[:,end] .= 0
        #left_edge =  ((σxx_minus_h2[:,1]-4*σxx_minus_h2[:,2]+3*σxx_minus_h2[:,3])./thick_edgeFace[:,1])/(2*dx)
        #σxx_minus_h2_x[:,end] = right_edge
        #σxx_minus_h2_x[:,1] = left_edge
        #σxx_minus_h2_x = hcat(zeros(grid.Ny+2,1),σxx_minus_h2_x,right_edge)
        #σxx_minus_h2_x = hcat(zeros(grid.Ny+2,1),σxx_minus_h2_x,zeros(grid.Ny+2,1))

        
        

        σyy_minus_h2_y = (σyy_minus_h2[3:end,:]-σyy_minus_h2[1:end-2,:])./thick_edgeFace[2:end-1,:]/(2*dx)
        σyy_minus_h2_y[1,:] .= (σyy_minus_h2[3,:]- σyy_minus_h2[2,:])/dx./thick_edgeFace[2,:]
        σyy_minus_h2_y[end,:] .= (σyy_minus_h2[end,:]- σyy_minus_h2[end-1,:])/dx./thick_edgeFace[end-1,:]#

        σxx_minus_h2_y = (σxx_minus_h2[3:end,:]-σxx_minus_h2[1:end-2,:])./thick_edgeFace[2:end-1,:]/(2*dx)
        σxx_minus_h2_y[1,:] .= (σxx_minus_h2[3,:]- σxx_minus_h2[2,:])/dx./thick_edgeFace[2,:]
        σxx_minus_h2_y[end,:] .= (σxx_minus_h2[end,:]- σxx_minus_h2[end-1,:])/dx./thick_edgeFace[end-1,:]
        #σyy_minus_h2_y[1,:] .= 0.0
        #σyy_minus_h2_y[end,:] .= 0.0

        #println(size(σxx_minus_h2_x[2:end-1,2:end-1]))
        #println(size(u[2:end-1,2:end-1]))
        #t = 0.5*dt^2*argand_upon_mach2.*(σxx[2:end-1,3:end]./η[2:end-1,3:end] - σxx[2:end-1,1:end-2]./η[2:end-1,1:end-2])/(2*dx)./thick_edgeFace[2:end-1,2:end-1]
        #####u[2:end-1,2:end-1] = u[2:end-1,2:end-1] +
            # dt*(σxx-h^2)_x
            #####dt/froude_number_squared.*σxx_minus_h2_x[2:end-1,:] +
            # dt*σxy,y
            ###dt/froude_number_squared./thick_edgeFace[2:end-1,2:end-1].*(σxy[3:end,2:end-1]-σxy[1:end-2,2:end-1])/(2*dx) +
            # 0.5*dt^2*(hu_x)_x
            #####0.5*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*hu_xx[2:end-1,2:end-1] + 
            # 0.5*dt^2*0.5*(hv_y)_x  
            ###0.25*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hv_y[:,3:end]-hv_y[:,1:end-2])/(2*dx) +
            #-0.5*dt^2*Ar*(σxx/η)_x
            #####-0.5*dt^2*argand_upon_mach2.*(σxx[2:end-1,3:end]./η[2:end-1,3:end] - σxx[2:end-1,1:end-2]./η[2:end-1,1:end-2])/(2*dx)./thick_edgeFace[2:end-1,2:end-1] ###+
            # 0.5*dt^2*(h*u_y)_y # Check if this shouldn't be thick_sy????
            ###0.125*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*hu_yy[2:end-1,2:end-1] +   
            # 0.5*dt^2*(h*v_x)_y
            ###0.125*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hv_x[3:end,:]-hv_x[1:end-2,:])/(2*dx) +
            # 0.5*dt^2*2*Ar*(σxy/η)_y
            ###-0.25*(dt^2)*argand_upon_mach2./thick_edgeFace[2:end-1,2:end-1].*(σxy[3:end,2:end-1]./η[3:end,2:end-1]-σxy[1:end-2,2:end-1]./η[1:end-2,2:end-1])/(2*dx)
        speed = sqrt.(u.^2+v.^2)[2:end-1,2:end-1]
        
        u[2:end-1,2:end-1] = u[2:end-1,2:end-1] + 
            # dt*(Fr*h)^(-1)*(σxx-h^2)_x, Factor of 2 from central difference
            0.5*(dt/dx)./(thick_edgeFace[2:end-1,2:end-1])./(froude_number_squared).*(σxx_minus_h2[2:end-1,3:end]-σxx_minus_h2[2:end-1,1:end-2]) +
            # 0.5*dt*(Fr*h)^(-1)*σxy,y, Factor of 2 from central difference
            0.25*(dt/dx)./thick_edgeFace[2:end-1,2:end-1]./froude_number_squared.*(σxy[3:end,2:end-1]-σxy[1:end-2,2:end-1]) +
            #-0.5*dt^2*(M^2*h)^(-1)*Ar*(σxx/η)_x, Factor of 2 from central difference
            -0.25*(dt^2/dx)./(thick_edgeFace[2:end-1,2:end-1]) * argand_upon_mach2.*(σxx[2:end-1,3:end]./η[2:end-1,3:end]-σxx[2:end-1,1:end-2]./η[2:end-1,1:end-2]) + 
            # 0.5*dt^2*(M^2*h)^(-1)*(hu_x)_x, using staggered second difference
            0.5*(dt/dx)^2/mach_number_squared.*(thick_sx[:,1:end-1].*u[2:end-1,1:end-2] - (thick_sx[:,1:end-1] + thick_sx[:,2:end]).*u[2:end-1,2:end-1] + thick_sx[:,2:end].*u[2:end-1,3:end])./thick_edgeFace[2:end-1,2:end-1]  +
            # 0.25*dt^2**(M^2*h)^(-1)*Ar*(σxy/η)_y, Factor of 2 from central differences
            -0.125*(dt^2/dx)*argand_upon_mach2./thick_edgeFace[2:end-1,2:end-1].*(σxy[3:end,2:end-1]./η[3:end,2:end-1]-σxy[1:end-2,2:end-1]./η[1:end-2,2:end-1]) +
            # 0.25*dt^2*(M^2*h)^(-1)*(h*v_y)_x 
            0.25*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hv_y[:,3:end]-hv_y[:,1:end-2])/(2*dx) +
            # 0.125*dt^2*(M^2*h)^(-1)*(h*v_x)_y 
            0.125*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hv_x[3:end,:]-hv_x[1:end-2,:])/(2*dx) +
            # 0.125*dt^2*(M^2*h)^(-1)*(h*u_y)_y
            0.125*(dt/dx)^2/mach_number_squared.*(thick_sy[1:end-1,:].*u[1:end-2,2:end-1] - (thick_sy[1:end-1,:] + thick_sy[2:end,:]).*u[2:end-1,2:end-1] + thick_sy[2:end,:].*u[3:end,2:end-1])./thick_edgeFace[2:end-1,2:end-1] +
            -dt*beta.*speed.^(1/3-1).*u[2:end-1,2:end-1]
        
            #unew = 0.25*(dt/dx)./thick_edgeFace[2:end-1,2:end-1]./froude_number_squared.*(σxy[3:end,2:end-1]-σxy[1:end-2,2:end-1]) 
        #println("Edge vel",unew[1,1])
        u[1,:] = u[2,:]
        u[end,:]=u[end-1,:]
        u[:,end] = u[:,end-1]
        u[:,1] = 2*U_inflow .- u[:,2] 
        #u[:,1]= 0 .- u[:,2]       # Average across cells = margin vel
        #u[:,end]= 0 .-u[:,end-1] # Average across cells = margin vel
        
        
        v[2:end-1,2:end-1] = v[2:end-1,2:end-1] + 
            # dt*(Fr*h)^(-1)*(σyy-h^2)_y, Factor of 2 from central difference
            # switched+
            0.5*(dt/dx)./(thick_edgeFace[2:end-1,2:end-1])./(froude_number_squared).*(σyy_minus_h2[3:end,2:end-1]-σyy_minus_h2[1:end-2,2:end-1]) +
            # 0.5*dt*(Fr*h)^(-1)*σxy,x, Factor of 2 from central difference
            0.25*(dt/dx)/froude_number_squared./thick_edgeFace[2:end-1,2:end-1].*(σxy[2:end-1,3:end]-σxy[2:end-1,1:end-2]) +
            #-0.5*dt^2*(M^2*h)^(-1)*Ar*(σyy/η)_y, Factor of 2 from central difference
            -0.25*(dt^2/dx)./(thick_edgeFace[2:end-1,2:end-1]) * argand_upon_mach2.*(σyy[3:end,2:end-1]./η[2:end-1,3:end]-σyy[1:end-2,2:end-1]./η[2:end-1,1:end-2]) + 
            # 0.5*dt^2*(M^2*h)^(-1)*(hv_y)_y, using staggered second difference
            0.5*(dt/dx)^2/mach_number_squared.*(thick_sy[1:end-1,:].*v[1:end-2,2:end-1] - (thick_sy[1:end-1,:] + thick_sy[2:end,:]).*v[2:end-1,2:end-1] + thick_sy[2:end,:].*v[3:end,2:end-1])./thick_edgeFace[2:end-1,2:end-1]  +
            # 0.5*dt^2**(M^2*h)^(-1)*Ar*(σxy/η)_x, Factor of 2 from central differences
            -0.125*(dt^2/dx)*argand_upon_mach2./thick_edgeFace[2:end-1,2:end-1].*(σxy[2:end-1,3:end]./η[2:end-1,3:end]-σxy[2:end-1,1:end-2]./η[2:end-1,1:end-2]) +
            # 0.25*dt^2*(M^2*h)^(-1)*(hu_x)_y 
            0.25*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hu_x[3:end,:]-hu_x[1:end-2,:])/(2*dx) +
            # 0.25*dt^2*(M^2*h)^(-1)*(h*u_y)_x 
            0.125*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hu_y[:,3:end,]-hu_y[:,1:end-2])/(2*dx) +
            # 0.25*dt^2*(M^2*h)^(-1)*(h*v_x)_x 
            0.125*(dt/dx)^2/mach_number_squared.*(thick_sx[:,1:end-1].*v[2:end-1,1:end-2] - (thick_sx[:,1:end-1] + thick_sx[:,2:end]).*v[2:end-1,2:end-1] + thick_sx[:,2:end].*v[2:end-1,3:end])./thick_edgeFace[2:end-1,2:end-1] +
            -dt*beta.*speed.^(1/3-1).*v[2:end-1,2:end-1]

        v[1,:] = - v[2,:].*mask .+ (1.0 .- mask).*v[2,:]
        v[end,:]=- v[end-1,:].*mask .+ (1.0 .- mask).*v[end-1,:]
        v[:,1]= v[:,2]       # Average across cells = margin vel
        v[:,end]= v[:,end-1] # Average across cells = margin vel
        
        

        
        σxx[2:end-1,2:end-1] = σxx[2:end-1,2:end-1] + 
            #(Fr/M)^2*(dt-dt^2*(F/M)^2*Ar/η*σxx)*(hu_x +0.5*hv_y -Ar/η*σxx)
            froude_upon_mach2*(dt .- 0.5*dt^2*froude_upon_mach2*argand_number./η[2:end-1,2:end-1]).*(hu_x[2:end-1,:] .+ 0.5*hv_y[:,2:end-1]  - argand_number*σxx[2:end-1,2:end-1]./η[2:end-1,2:end-1]) +
            # 0.5*dt^2/M^2*h*[(1/h)*(σxx -h^2)_x]_x, staggered second difference
            + 0.5*(dt/dx)^2/mach_number_squared*(σxx_minus_h2[2:end-1,1:end-2]./thick_sx[:,1:end-1] .- (1 ./thick_sx[:,1:end-1] .+ 1 ./thick_sx[:,2:end]).*σxx_minus_h2[2:end-1,2:end-1] + σxx_minus_h2[2:end-1,3:end]./thick_sx[:,2:end]).*thick_edgeFace[2:end-1,2:end-1] +
            # 0.25*dt^2/M^2*h*[(1/h)*(σyy -h^2)_y]_y, staggered second difference
            + 0.25*(dt/dx)^2/mach_number_squared*(σyy_minus_h2[1:end-2,2:end-1]./thick_sy[1:end-1,:] .- (1 ./thick_sy[1:end-1,:] .+ 1 ./thick_sy[2:end,:]).*σyy_minus_h2[2:end-1,2:end-1] + σyy_minus_h2[3:end,2:end-1]./thick_sy[2:end,:]).*thick_edgeFace[2:end-1,2:end-1] +
            # 0.25*dt^2/M^2*h*[(1/h)*σxy_y]_x,
            + 0.25*dt^2/mach_number_squared.*thick_edgeFace[2:end-1,2:end-1].*(σxy_y[:,3:end]-σxy_y[:,1:end-2])/(2*dx) +
            # 0.25*dt^2/M^2*h*[(1/h)*σxy_x]_y,
            + 0.125*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy_x[3:end,:]-σxy_x[1:end-2,:])/(2*dx)
       
        σxx[1,:]   = σxx[2,:]
        σxx[end,:] =σxx[end-1,:]
        σxx[:,1] = σxx[:,2]
        σxx[:,end] = 2.0.*(0.5.*(thick_edgeFace[:,end]+thick_edgeFace[:,end-1])).^2  .- σxx[:,end-1]        
        

            
        σyy[2:end-1,2:end-1] = σyy[2:end-1,2:end-1] + 
            #(Fr/M)^2*(dt-dt^2*(F/M)^2*Ar/η*σyy)*(hv_y +0.5*hu_x -Ar/η*σyy)
            froude_upon_mach2*(dt .- 0.5*dt^2*froude_upon_mach2*argand_number./η[2:end-1,2:end-1]).*(hv_y[:,2:end-1] .+ 0.5*hu_x[2:end-1,:]  - argand_number*σyy[2:end-1,2:end-1]./η[2:end-1,2:end-1]) +
            # 0.5*dt^2/M^2*h*[(1/h)*(σyy -0.5*h^2)_y]_y, staggered second difference
            + 0.5*(dt/dx)^2/mach_number_squared*(σyy_minus_h2[1:end-2,2:end-1]./thick_sy[1:end-1,:] .- (1 ./thick_sy[1:end-1,:] .+ 1 ./thick_sy[2:end,:]).*σyy_minus_h2[2:end-1,2:end-1] + σyy_minus_h2[3:end,2:end-1]./thick_sy[2:end,:]).*thick_edgeFace[2:end-1,2:end-1] +
            # 0.25*dt^2/M^2*h*[(1/h)*(σxx -0.5*h^2)_x]_x, staggered second difference
            + 0.25*(dt/dx)^2/mach_number_squared*(σxx_minus_h2[2:end-1,1:end-2]./thick_sx[:,1:end-1] .- (1 ./thick_sx[:,1:end-1] .+ 1 ./thick_sx[:,2:end]).*σxx_minus_h2[2:end-1,2:end-1] + σxx_minus_h2[2:end-1,3:end]./thick_sx[:,2:end]).*thick_edgeFace[2:end-1,2:end-1] +
            # 0.25*dt^2/M^2*h*[(1/h)*σxy_x]_y,
            + 0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy_x[3:end,:]-σxy_x[1:end-2,:])/(2*dx) +
            # 0.125*dt^2/M^2*h*[(1/h)*σxy_y]_x,
            + 0.125*dt^2/mach_number_squared.*thick_edgeFace[2:end-1,2:end-1].*(σxy_y[:,3:end]-σxy_y[:,1:end-2])/(2*dx) 
        
        σyy[1,:] =  σyy[2,:].*mask .+ (1.0 .- mask).*(2.0.*(0.5.*(thick_edgeFace[1,:]+thick_edgeFace[2,:])).^2  .- σyy[2,:])    
        σyy[end,:] =  σyy[end-1,:].*mask .+ (1.0 .- mask).*(2.0.*(0.5.*(thick_edgeFace[end,:]+thick_edgeFace[end-1,:])).^2  .- σyy[end-1,:])
       
        σyy[:,end]=  σyy[:,end-1]
        σyy[:,1]= σyy[:,2]
        
        
        h_mid_x = 0.5*(thick_edgeFace[:,1:end-1]+thick_edgeFace[:,2:end])
        h_mid_y = 0.5*(thick_edgeFace[1:end-1,:]+thick_edgeFace[2:end,:])
      
        σxy[2:end-1,2:end-1] = σxy[2:end-1,2:end-1] + 
            #dt*(Fr/M)^2*(0.5*dt-0.25*dt*(F/M)^2*Ar/η)*(0.5*(hu_y + hv_x) -Ar*σxy/eta)
            froude_upon_mach2*(dt .- 0.5*dt^2*froude_upon_mach2*argand_number./η[2:end-1,2:end-1]).*(0.5*(hu_y[:,2:end-1] + hv_x[2:end-1,:]) - argand_number*σxy[2:end-1,2:end-1]./η[2:end-1,2:end-1]) +
            #0.25*dt^2/M^2*h*[(1/h)(σxx-h^2)_x]_y
            0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxx_minus_h2_x[3:end,:]-σxx_minus_h2_x[1:end-2,:])/(2*dx) +
            #0.25*dt^2/M^2*h*[(1/h)(σyy-h^2)_y]_x
            0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σyy_minus_h2_y[:,3:end]-σyy_minus_h2_y[:,1:end-2])/(2*dx) +
            #0.125*dt^2/M^2*h*[(1/h)σxy_y]_y ## Why is there an extra factor of 2 in the ice thickness average in the 2nd derivative????????
            0.125*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy[2:end-1,1:end-2]./thick_sx[:,1:end-1] .- (1 ./thick_sx[:,1:end-1] .+ 1 ./thick_sx[:,2:end]).*σxy[2:end-1,2:end-1] + σxy[2:end-1,3:end]./thick_sx[:,2:end])/(dx^2) +
            #0.125*0.125*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy[2:end-1,1:end-2]./h_mid_x[2:end-1,1:end-1] -(1.0./h_mid_x[2:end-1,1:end-1] + 1.0./h_mid_x[2:end-1,2:end]).*σxy[2:end-1,2:end-1] + σxy[2:end-1,3:end]./h_mid_x[2:end-1,2:end])/(dx^2) +
            #0.125*dt^2/M^2*h*[(1/h)σxy_x]_x
            0.125*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy[1:end-2,2:end-1]./thick_sy[1:end-1,:] .- (1 ./thick_sy[1:end-1,:] .+ 1 ./thick_sy[2:end,:]).*σxy[2:end-1,2:end-1] + σxy[3:end,2:end-1]./thick_sy[2:end,:])/(dx^2)
            #0.125*0.125*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy[1:end-2,2:end-1]./h_mid_y[1:end-1,2:end-1] -(1.0./h_mid_y[1:end-1,2:end-1] + 1.0./h_mid_y[2:end,2:end-1]).*σxy[2:end-1,2:end-1] + σxy[3:end,2:end-1]./h_mid_y[2:end,2:end-1])/(dx^2)
        σxy[:,end] = -σxy[:,end-1]
        σxy[2,:] = -σxy[end-1,:]
        #tau_m1 = 1.75*(1.1*u[1,:]).*thick_edgeFace[1,:]
        #au_m2 = 1.75*(1.1*u[end,:]).*thick_edgeFace[end,:]
        tau_m1 = tau_margin.*thick_edgeFace[1,:]
        tau_m2 = tau_margin.*thick_edgeFace[end,:]
        σxy[1,:] =   2*tau_m1.*mask .- σxy[2,:]
        σxy[end,:]= -2*tau_m1.*mask .- σxy[end-1,:]
        

        #σxy[:,1] =  2*tau_margin .- σxy[:,2]
        #σxy[:,end]= 2*tau_margin .- σxy[:,end-1]

        
        
        T.xx[1:end,1:end] =  (σxx-0.5*σyy)/0.75
        T.yy[1:end,1:end] =  (σyy-0.5*σxx)/0.75
        T.xy[1:end,1:end] =   σxy
        #T.xx[:,end]=T.xx[:,end-1]+(T.xx[:,end-3]-4*T.xx[:,end-2]+3*T.xx[:,end-1])
        #T.yy[:,end]=T.yy[:,end-1]+(T.yy[:,end-3]-4*T.yy[:,end-2]+3*T.yy[:,end-1])       
        t += dt  
             
    end
    #vel.u = copy(u)
    #vel.v = copy(v)
    return t
end


#function update_thick(vel,T,h,U_inflow,V_inflow,H_inflow,accum,dx,dt,dim;type="superbee")
#    #h_old = copy(h)
#    #advect_LW(vel,h,accum,dx,dt,type=type)
#    h[:,1] .= 2*H_inflow .- h[:,2]
#    h[:,end] .= h[:,end-1]
#    h[end,:] .= h[end-1,:]
#    h[1,:] .= h[2,:]
#    h[h.<0.001] .= 0.001
#    return 
#end
#=
function update(vel,T,h,U_inflow,V_inflow,H_inflow,accum,grid,dim;num_steps=1,num_it=1,mu=1.0)
    t=0
    for timeCount in 1:num_steps 
        dt=update_vel(vel,T,h,U_inflow,V_inflow,grid,dim,num_steps=num_it,mu=mu)
        # This doesn't seem to work :(
        #=
        h_new=advect_staggered_grid(vel,h,accum,grid.dx,grid.dt)
        if length(H_inflow)>1
            H_inflow_staggered = 0.5*(H_inflow[1:end-1]+H_inflow[2:end])
        else
            H_inflow_staggered = H_inflow
        end

        # Set boundary conditions on staggered grid
        h_new[:,1] .= H_inflow_staggered
        h_new[:,end] = h_new[:,end-1]
        h_new[1,:] = h_new[2,:]
        h_new[end,:] = h_new[end-1,:]
        
        # Interpolate back to regular grid
        h[2:end-1,2:end-1] = cellCenter(h_new)
        h[:,1]= 2*H_inflow .- h[:,2]
        h[:,end] .= h[:,end-1]
        h[end,:] .= h[end-1,:]
        h[:,end]= h[:,end-1] .+ (h[:,end-1]-h[:,end-2])
        h[1,:] .= h[2,:]
        h[h.<0.01] .= 0.01
        =#
        # This does seem to work :)
        update_thick(vel,T,h,U_inflow,V_inflow,H_inflow,accum,grid.dx,grid.dt,dim)
        h[:,end]= h[:,end-1] .+ (h[:,end-1]-h[:,end-2])
        t+=dt
    end
    return t
end
=#

function update_vel_LW(vel,T,h,U_inflow,V_inflow,grid,prms,visc_prms;num_steps=1,mu=1.0,type="superbee",drag_x=0.0,drag_y=0.0,tau_margin=0.0)
    u = vel.u 
    v = vel.v 

    txx = T.xx 
    tyy = T.yy 
    txy = T.xy

    
    dt = grid.dt 
    dx = grid.dx
    froude_number_squared = prms.froude_number_squared
    argand_upon_mach2 = prms.argand_upon_mach2
    mach_number_squared = prms.mach_number_squared
    froude_upon_mach2 = prms.froude_upon_mach2
    argand_number = prms.argand_number
    mach_number = prms.mach_number


    # Square and ratio of some non-dimensional parameters
    #froude_number_squared = froude_number^2
    #mach_number_squared = mach_number^2
    #froude_upon_mach2 = froude_number_squared/mach_number_squared
    #argand_upon_mach2 = argand_number / mach_number_squared


    thick_edgeFace = h
    thick = cellCenter(h)
    thick_cellCenter = thick
    #thick = ice_thick.thick_cellCenter
    #thick_grad_x,thick_grad_y = grad_thick(thick_cellCenter,grid)
    
    thick_sy =  0.5*(thick_cellCenter[:,2:end]+thick_cellCenter[:,1:end-1])
    thick_sx =  0.5*(thick_cellCenter[2:end,:]+thick_cellCenter[1:end-1,:])
    #thick_sx,thick_sy = cellCenter(ice_thick.thick_cellCenter)
    #h_end=0.5*(thick_cellCenter[1:end-1,end]+thick_cellCenter[2:end,end])

    set_bc(u,v,T.xx,T.yy,T.xy,thick,U_inflow,V_inflow)


    t=0
    σxx = T.xx + 0.5*T.yy
    σyy = T.yy + 0.5*T.xx
    σxy = 0.5*T.xy


    
    
    for timeCount in 1:num_steps
       
        η = mu.*visc(T,thick_edgeFace,visc_prms)

        
        
        
        
        # Update time
       
        #println(t)
        #σxy = 0.5*txy
        #set_bc(u,v,T.xx,T.yy,T.xy,thick,U_inflow,V_inflow)

        #σyy = T.yy + 0.5*T.xx
        #σxx = T.xx + 0.5*T.yy 
     
       
    
        
        # Calculate viscosity
        σxx_minus_h2 = σxx- thick_edgeFace.^2
        σyy_minus_h2 = σyy- thick_edgeFace.^2

        # Calculate diffusion terms (hu_x)_x, (hu_y)_y
        hu_xx,hu_yy = diffusivity(thick_edgeFace,u,dx)
        # Calculate diffusion terms (hv_x)_x, (hv_y)_y
        hv_xx,hv_yy = diffusivity(thick_edgeFace,v,dx)

       
         # Calculate [σyy_x)/thick]_x,[σyy_y)/thick]_y
         #σxx_minus_h2_xx,σxx_minus_h2_yy = diffusivity(1.0./thick_edgeFace,σxx_minus_h2,dx)
         # Calculate [σyy_x)/thick]_x,[σyy_y)/thick]_y
         #σyy_minus_h2_xx,σyy_minus_h2_yy = diffusivity(1.0./thick_edgeFace,σyy_minus_h2,dx)



        hv_y = (v[3:end,:]-v[1:end-2,:]).*thick_edgeFace[2:end-1,:]/(2*dx)
        hu_y = (u[3:end,:]-u[1:end-2,:]).*thick_edgeFace[2:end-1,:]/(2*dx)
        hu_x = (u[:,3:end]-u[:,1:end-2]).*thick_edgeFace[:,2:end-1]/(2*dx)
        hv_x = (v[:,3:end]-v[:,1:end-2]).*thick_edgeFace[:,2:end-1]/(2*dx)

        
        σxy_y   = (σxy[3:end,:]-σxy[1:end-2,:])./thick_edgeFace[2:end-1,:]/(2*dx) 
        σxy_x   = (σxy[:,3:end]-σxy[:,1:end-2])./thick_edgeFace[:,2:end-1]/(2*dx)

       

        #σxx_minus_h2_x = (σxx_minus_h2[:,2:end]-σxx_minus_h2[:,1:end-1])./(0.5*(thick_edgeFace[:,1:end-1]+thick_edgeFace[:,2:end]))
        #right_edge = (σxx_minus_h2[:,end-2]-4*σxx_minus_h2[:,end-1]+3*σxx_minus_h2[:,end])./thick_edgeFace[:,end]
        #σxx_minus_h2_x = hcat(zeros(grid.Ny+2,1),σxx_minus_h2_x,right_edge)/(dx)


        #σxx_minus_h2_x[:,1] .= 0.0 # Assume gradient vanishes at inflow boundary . . .
        #σxx_minus_h2_xx = (σxx_minus_h2_x[:,2:end]-σxx_minus_h2_x[:,1:end-1])/(dx)

        σxx_minus_h2_x = (σxx_minus_h2[:,3:end]-σxx_minus_h2[:,1:end-2])./thick_edgeFace[:,2:end-1]/(2*dx)

        #σxx_minus_h2_x = (σxx_minus_h2[:,2:end-1]-σxx_minus_h2[:,1:end-2])./thick_edgeFace[:,2:end-1]/(dx)

        σxx_minus_h2_x[:,1] .= (σxx_minus_h2[:,3]- σxx_minus_h2[:,2])/dx./thick_edgeFace[:,2]
        σxx_minus_h2_x[:,end] .= (σxx_minus_h2[:,end]- σxx_minus_h2[:,end-1])/dx./thick_edgeFace[:,end-1]
        #σxx_minus_h2_x[:,1] = σxx_minus_h2_x[:,2]
        
        #σxx_minus_h2_x[:,1] .= 0
        #right_edge = ((σxx_minus_h2[:,end-2]-4*σxx_minus_h2[:,end-1]+3*σxx_minus_h2[:,end])./thick_edgeFace[:,end-1])/(2*dx)
        #σxx_minus_h2_x[:,end] = right_edge
        #σxx_minus_h2_x[:,end] .= 0
        #left_edge =  ((σxx_minus_h2[:,1]-4*σxx_minus_h2[:,2]+3*σxx_minus_h2[:,3])./thick_edgeFace[:,1])/(2*dx)
        #σxx_minus_h2_x[:,end] = right_edge
        #σxx_minus_h2_x[:,1] = left_edge
        #σxx_minus_h2_x = hcat(zeros(grid.Ny+2,1),σxx_minus_h2_x,right_edge)
        #σxx_minus_h2_x = hcat(zeros(grid.Ny+2,1),σxx_minus_h2_x,zeros(grid.Ny+2,1))

        
        

        σyy_minus_h2_y = (σyy_minus_h2[3:end,:]-σyy_minus_h2[1:end-2,:])./thick_edgeFace[2:end-1,:]/(2*dx)
        σyy_minus_h2_y[1,:] .= (σyy_minus_h2[3,:]- σyy_minus_h2[2,:])/dx./thick_edgeFace[2,:]# Assume gradient vanishes at inflow boundary . . .
        σyy_minus_h2_y[end,:] .= (σyy_minus_h2[end,:]- σyy_minus_h2[end-1,:])/dx./thick_edgeFace[end-1,:]#
       

        w1 = 0.5*(+froude_number_squared.*thick_edgeFace.*u/mach_number+σxx)
        w2 = 0.5*(-froude_number_squared.*thick_edgeFace.*u/mach_number+σxx)
               
        source_x =  -0.5*froude_number_squared.*u[2:end-1,2:end-1].*(thick_edgeFace[2:end-1,3:end].-thick_edgeFace[2:end-1,1:end-2])/dx/2/mach_number_squared .+
                    -0.5*froude_number_squared.*thick_edgeFace[2:end-1,2:end-1].*(thick_edgeFace[2:end-1,3:end].-thick_edgeFace[2:end-1,1:end-2])/dx/2/mach_number .+
                    -0.5*argand_number.*σxx[2:end-1,2:end-1]./η[2:end-1,2:end-1]*froude_number_squared/mach_number_squared .+
                    -0.5*froude_number_squared.*drag_x/mach_number #+
                    #-0.5*argand_number*froude_upon_mach2*(1 .- 0.5*dt*froude_upon_mach2*argand_number./η[2:end-1,2:end-1]).*σxx[2:end-1,2:end-1]./η[2:end-1,2:end-1]
    
        source_y =  -0.5*froude_number_squared.*u[2:end-1,2:end-1].*(thick_edgeFace[2:end-1,3:end].-thick_edgeFace[2:end-1,1:end-2])/dx/2/mach_number_squared .+
                    +0.5*froude_number_squared.*thick_edgeFace[2:end-1,2:end-1].*(thick_edgeFace[2:end-1,3:end].-thick_edgeFace[2:end-1,1:end-2])/dx/2/mach_number .+
                    -0.5*argand_number.*σxx[2:end-1,2:end-1]./η[2:end-1,2:end-1]*froude_number_squared/mach_number_squared .+
                    +0.5*froude_number_squared.*drag_x/mach_number

        #source_x[:,end] .= 0.0
        #source_y[:,end] .= 0.0

        #source_x[:,1] .= source_x[:,2]
        #source_y[:,1] .= source_y[:,2]
           

        #c = -ones(size(u))/mach_number
        c1 = 0.0*u .- 1/mach_number
        vel_seismic1 = Velocity(c1,c1)
        Fx1,junk = flux(vel_seismic1,w1,dx,dt;type=type)
        c2 = 0.0*u .+ 1/mach_number
        vel_seismic2 = Velocity(c2,c2)
        Fx2,junk = flux(vel_seismic2,w2,dx,dt;type=type)
        #Fx1[:,end] .= Fx1[:,end-1]
        #Fx2[:,end] .= Fx2[:,end-1]
        #Fx2[:,2] .= Fx2[:,1]
        #Fx2[:,end] .= Fx2[:,end-1]
        

        
        w1[2:end-1,2:end-1] .= w1[2:end-1,2:end-1] .- dt/dx*(Fx1[:,2:end]-Fx1[:,1:end-1]) .+
                    source_x.*dt
                    #dt/dx.*(s1[:,2:end]-s1[:,1:end-1])
        w2[2:end-1,2:end-1] .= w2[2:end-1,2:end-1] .- dt/dx*(Fx2[:,2:end]-Fx2[:,1:end-1]) .+
                    source_y.*dt 
                    #dt/dx.*(s2[:,2:end]-s2[:,1:end-1])

        #w1[end,:] .=w1[end-1,:]# - w1[end-2,:]
        #w2[end,:] .= w2[end-1,:]# - w2[end-2,:]
        #w1[1,:] .= 2*w1[2,:] - w1[2,:]
        #w2[1,:] .= 2*w2[2,:] - w2[2,:]
        #w1[:,1] = w1[:,2]
        #w2[:,1]= w2[:,2]
        
        u[:,:] .=  mach_number./froude_number_squared./thick_edgeFace.*(w1 .- w2)
        σxx[:,:] .= w1 .+ w2

        #u[1,:] = u[2,:]
        u[end,:]=u[end-1,:]
        u[:,end] = u[:,end-1]
        u[:,1]= 2*U_inflow .- u[:,2] 

        σxx[1,:]  = σxx[2,:]
        σxx[end,:] =σxx[end-1,:]
        σxx[2:end,end] = 2.0.*thick[:,end].^2  .- σxx[2:end,end-1]        
        σxx[:,1] = σxx[:,2]
        
        #=
        v[2:end-1,2:end-1] = v[2:end-1,2:end-1] + 
            dt/froude_number_squared.*σyy_minus_h2_y[:,2:end-1] +
            dt/froude_number_squared./thick_edgeFace[2:end-1,2:end-1].*(σxy[2:end-1,3:end]-σxy[2:end-1,1:end-2])/(2*dx) +
            0.5*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*hv_yy[2:end-1,2:end-1] +        
            0.25*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hu_x[3:end,:]-hu_x[1:end-2,:])/(2*dx) -
            0.5*dt^2*argand_upon_mach2.*(σyy[3:end,2:end-1]./η[3:end,2:end-1] - σyy[1:end-2,2:end-1]./η[1:end-2,2:end-1])/(2*dx)./thick_edgeFace[2:end-1,2:end-1] +            
            0.125*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*hv_xx[2:end-1,2:end-1] +   
            0.125*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hu_y[:,3:end]-hu_y[:,1:end-2])/(2*dx) -
            0.25*dt^2*argand_upon_mach2./thick_edgeFace[2:end-1,2:end-1].*(σxy[2:end-1,3:end]./η[2:end-1,3:end]-σxy[2:end-1,1:end-2]./η[2:end-1,1:end-2])/(2*dx)
    
        v[1,:] = 2*V_inflow .- v[2,:]
        v[end,:]=2*V_inflow .- v[end-1,:]
        v[:,1]= v[:,2]       # Average across cells = margin vel
        v[:,end]= v[:,end-1] # Average across cells = margin vel

        σyy[2:end-1,2:end-1] = σyy[2:end-1,2:end-1] + 
            (dt*froude_upon_mach2 .- 0.5*dt^2*froude_upon_mach2^2*argand_number./η[2:end-1,2:end-1]).*(hv_y[:,2:end-1] .+0.5*hu_x[2:end-1,:]  - argand_number*σyy[2:end-1,2:end-1]./η[2:end-1,2:end-1]) +
            0.5*dt^2/mach_number_squared*(σyy_minus_h2[1:end-2,2:end-1]./thick_sy[1:end-1,:] .- (1 ./thick_sy[1:end-1,:] .+ 1 ./thick_sy[2:end,:]).*σyy_minus_h2[2:end-1,2:end-1] + σyy_minus_h2[3:end,2:end-1]./thick_sy[2:end,:])/(dx^2).*thick_edgeFace[2:end-1,2:end-1] +    
            #0.5*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*σyy_minus_h2_yy[2:end-1,2:end-1] +    
            0.5*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy_x[3:end,:]-σxy_x[1:end-2,:])/(2*dx) +
            0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxx_minus_h2[2:end-1,1:end-2]./thick_sx[:,1:end-1] .- (1 ./thick_sx[:,1:end-1] .+ 1 ./thick_sx[:,2:end]).*σxx_minus_h2[2:end-1,2:end-1] + σxx_minus_h2[2:end-1,3:end]./thick_sx[:,2:end])/dx^2 +
            #0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*σxx_minus_h2_xx[2:end-1,2:end-1] +
            0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy_y[:,3:end]-σxy_y[:,1:end-2])/(2*dx)
        
        σyy[1,:] =  σyy[2,:]
        σyy[end,:] =  σyy[end-1,:]
        σyy[:,end]=  σyy[:,end-1]
        σyy[:,1]= σyy[:,2]
            # Update the rest of the velocity with shear terms
        #u[2:end-1,2:end-1] = u[2:end-1,2:end-1] .+
         #   dt/froude_number_squared./thick_edgeFace[2:end-1,2:end-1].*(σxy[3:end,2:end-1]-σxy[1:end-2,2:end-1])/(2*dx) #+
            #0.25*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hv_y[:,3:end]-hv_y[:,1:end-2])/(2*dx) +
            #0.125*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*hu_yy[2:end-1,2:end-1] +   
            #0.125*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hv_x[3:end,:]-hv_x[1:end-2,:])/(2*dx) +
            #-0.25*(dt^2)*argand_upon_mach2./thick_edgeFace[2:end-1,2:end-1].*(σxy[3:end,2:end-1]./η[3:end,2:end-1]-σxy[1:end-2,2:end-1]./η[1:end-2,2:end-1])/(2*dx)


        

        u[1,:] = u[2,:]
        u[end,:]=u[end-1,:]
        u[:,end] = u[:,end-1]
        u[:,1]=2*U_inflow .- u[:,2] 
        
        

        #σxx[2:end-1,2:end-1] = σxx[2:end-1,2:end-1] + 
         #   (dt*froude_upon_mach2 .- 0.5*dt^2*froude_upon_mach2^2*argand_number./η[2:end-1,2:end-1]).*(0.5*hv_y[:,2:end-1]) #+
            #0.5*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy_y[:,3:end]-σxy_y[:,1:end-2])/(2*dx) +
            #0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy_x[3:end,:]-σxy_x[1:end-2,:])/(2*dx)

        σxx[1,:]  = σxx[2,:]
        σxx[end,:] =σxx[end-1,:]
        σxx[2:end,end] = 2.0.*thick[:,end].^2  .- σxx[2:end,end-1]        
        σxx[:,1] = σxx[:,2]
        =#
        
        
        
        w1 = 0.5*(froude_number_squared.*thick_edgeFace.*v/mach_number+σyy)
        w2 = 0.5*(-froude_number_squared.*thick_edgeFace.*v/mach_number+σyy)
        
        source_x =  -0.5*froude_number_squared.*v[2:end-1,2:end-1].*(thick_edgeFace[3:end,2:end-1].-thick_edgeFace[1:end-2,2:end-1])/dx/2/mach_number_squared .+
                    -0.5*froude_number_squared.*thick_edgeFace[2:end-1,2:end-1].*(thick_edgeFace[3:end,2:end-1].-thick_edgeFace[1:end-2,2:end-1])/dx/2/mach_number .+
                    -0.5*argand_number.*σyy[2:end-1,2:end-1]./η[2:end-1,2:end-1]*froude_number_squared/mach_number_squared .+
                    +0.5*froude_number_squared.*thick_edgeFace[2:end-1,2:end-1].*drag_y/mach_number

        source_y =  -0.5*froude_number_squared.*v[2:end-1,2:end-1].*(thick_edgeFace[3:end,2:end-1].-thick_edgeFace[1:end-2,2:end-1])/dx/2/mach_number_squared .+
                    +0.5*froude_number_squared.*thick_edgeFace[2:end-1,2:end-1].*(thick_edgeFace[3:end,2:end-1].-thick_edgeFace[1:end-2,2:end-1])/dx/2/mach_number .+
                    -0.5*argand_number.*σyy[2:end-1,2:end-1]./η[2:end-1,2:end-1]*froude_number_squared/mach_number_squared .+
                    +0.5*froude_number_squared.*thick_edgeFace[2:end-1,2:end-1].*drag_y/mach_number
        #source_x[end,:] .= 0.0
        #source_y[end,:] .= 0.0
        #source_x[1,:] .= 0.0
        #source_y[1,:] .= 0.0
           
        
        #c = -ones(size(u))/mach_number
        c1 = 0.0*v .- 1/mach_number
        vel_seismic1 = Velocity(c1,c1)
        junk,Fx1 = flux(vel_seismic1,w1,dx,dt;type=type)
        c2 = 0.0*v .+ 1/mach_number
        vel_seismic2 = Velocity(c2,c2)
        junk,Fx2 = flux(vel_seismic2,w2,dx,dt;type=type)
        #Fx1[end,:] .= Fx1[end-1,:]  
        #Fx1[1,:] .= Fx1[2,:] 

        w1[2:end-1,2:end-1] .= w1[2:end-1,2:end-1] .- dt/dx*(Fx1[2:end,:]-Fx1[1:end-1,:])# .+
                    source_x.*dt
                    #dt/dx.*(s1[:,2:end]-s1[:,1:end-1])

        w2[2:end-1,2:end-1] .= w2[2:end-1,2:end-1] .- dt/dx*(Fx2[2:end,:]-Fx2[1:end-1,:])# .+
                    source_y.*dt 

        v[:,:] .=  mach_number./froude_number_squared./thick_edgeFace.*(w1 .- w2)
        σyy[:,:] .= w1 .+ w2


        v[1,:] = 2*V_inflow .- v[2,:]
        v[end,:]=2*V_inflow .- v[end-1,:]
        v[:,1]= v[:,2]       # Average across cells = margin vel
        v[:,end]= v[:,end-1] # Average across cells = margin vel

        σyy[1,:] =  σyy[2,:]
        σyy[end,:] =  σyy[end-1,:]
        σyy[:,end]=  σyy[:,end-1]
        σyy[:,1]= σyy[:,2]

        #=
        
       
        #v[2:end-1,2:end-1] = v[2:end-1,2:end-1] + 
        #    dt/froude_number_squared./thick_edgeFace[2:end-1,2:end-1].*(σxy[2:end-1,3:end]-σxy[2:end-1,1:end-2])/(2*dx) #+
            #0.25*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hu_x[3:end,:]-hu_x[1:end-2,:])/(2*dx) -
            #0.125*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*hv_xx[2:end-1,2:end-1] +   
            #0.125*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hu_y[:,3:end]-hu_y[:,1:end-2])/(2*dx) -
            #0.25*dt^2*argand_upon_mach2./thick_edgeFace[2:end-1,2:end-1].*(σxy[2:end-1,3:end]./η[2:end-1,3:end]-σxy[2:end-1,1:end-2]./η[2:end-1,1:end-2])/(2*dx)
    
       

        #v[1,:] = 2*V_inflow .- v[2,:]
        ##v[end,:]=2*V_inflow .- v[end-1,:]
        #v[:,1]= v[:,2]       # Average across cells = margin vel
        #v[:,end]= v[:,end-1] # Average across cells = margin vel

        
        
        
      
    
        
        #σyy[2:end-1,2:end-1] = σyy[2:end-1,2:end-1] + 
        #    (dt*froude_upon_mach2 .- 0.5*dt^2*froude_upon_mach2^2*argand_number./η[2:end-1,2:end-1]).*(0.5*hu_x[2:end-1,:]) #+
            #0.5*dt^2/mach_number_squared*(σyy_minus_h2[1:end-2,2:end-1]./thick_sy[1:end-1,:] .- (1 ./thick_sy[1:end-1,:] .+ 1 ./thick_sy[2:end,:]).*σyy_minus_h2[2:end-1,2:end-1] + σyy_minus_h2[3:end,2:end-1]./thick_sy[2:end,:])/(dx^2).*thick_edgeFace[2:end-1,2:end-1] +    
            #0.5*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy_x[3:end,:]-σxy_x[1:end-2,:])/(2*dx) +
            #0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxx_minus_h2[2:end-1,1:end-2]./thick_sx[:,1:end-1] .- (1 ./thick_sx[:,1:end-1] .+ 1 ./thick_sx[:,2:end]).*σxx_minus_h2[2:end-1,2:end-1] + σxx_minus_h2[2:end-1,3:end]./thick_sx[:,2:end])/dx^2 +
            #0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy_y[:,3:end]-σxy_y[:,1:end-2])/(2*dx)
        
        #σyy[1,:] =  σyy[2,:]
        #σyy[end,:] =  σyy[end-1,:]
        #σyy[:,end]=  σyy[:,end-1]
        #σyy[:,1]= σyy[:,2]

        
        #σyy[:,end]= h_end[1].^2 .- σyy[:,end-1] 
        =#
        
        
        # Double check size of term????      
        σxy[2:end-1,2:end-1] = σxy[2:end-1,2:end-1] + 
            (0.5*dt*froude_upon_mach2 .- 0.25*dt^2*froude_upon_mach2^2*argand_number./η[2:end-1,2:end-1]).*(0.5*(hu_y[:,2:end-1] + hv_x[2:end-1,:]) - 2*argand_number*σxy[2:end-1,2:end-1]./η[2:end-1,2:end-1]) +
            0.125*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxx_minus_h2_x[3:end,:]-σxx_minus_h2_x[1:end-2,:])/(2*dx) +
            0.125*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σyy_minus_h2_y[:,3:end]-σyy_minus_h2_y[:,1:end-2])/(2*dx) +
             0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy[2:end-1,1:end-2]./thick_sx[:,1:end-1] .- (1 ./thick_sx[:,1:end-1] .+ 1 ./thick_sx[:,2:end]).*σxy[2:end-1,2:end-1] + σxy[2:end-1,3:end]./thick_sx[:,2:end])/(dx^2) +
             0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy[1:end-2,2:end-1]./thick_sy[1:end-1,:] .- (1 ./thick_sy[1:end-1,:] .+ 1 ./thick_sy[2:end,:]).*σxy[2:end-1,2:end-1] + σxy[3:end,2:end-1]./thick_sy[2:end,:])/(dx^2)
        σxy[:,end] = -σxy[:,end-1];σxy[:,1]  = σxy[:,2]
        #σxy[1,:]   = -σxy[2,:];    σxy[end,:]= -σxy[end-1,:] 

        mask = ones(grid.Nx)
        σxy[1,:] =   2*tau_margin .- σxy[2,:]
        σxy[end,:]= -2*tau_margin .- σxy[end-1,:]

        #σxy[1,:] =   tau_margin 
        #σxy[end,:]= -tau_margin 
            
       
        #σxy=σxy/2
       
        T.xx[1:end,1:end] =  (σxx-0.5*σyy)/0.75
        T.yy[1:end,1:end] =  (σyy-0.5*σxx)/0.75
        T.xy[1:end,1:end] = σxy
        #σyy = T.yy + 0.5*T.xx

        #T.yy[1,:] =  T.yy[2,:]
        #T.yy[end,:]= T.yy[end-1,:]
        #T.yy[:,end]= T.yy[:,end-1]
        #T.yy[:,1]= T.yy[:,2]
        #T.yy = T.yy .- Statistics.mean(T.yy)
        
        #σyy = T.yy + 0.5*T.xx
        #σxx = T.xx + 0.5*T.yy
        
        
        t += dt  
             
    end
    return t
end


function update_thick(vel,T,h,U_inflow,V_inflow,H_inflow,accum,dx,dt,dim;type="superbee")
    #h_old = copy(h)
    advect_LW(vel,h,accum,dx,dt,type=type)
    #advect(vel,h,accum,dx,dt)
    h[end,:] .= h[end-1,:]
    h[1,:] .= h[2,:]
    h[:,1] .= 2*H_inflow .- h[:,2]
    h[:,end] .= h[:,end-1]
   
    h[h.<0.001] .= 0.001
    return 
end

function update(vel,T,h,U_inflow,V_inflow,H_inflow,accum,grid,dim,visc_prms;num_steps=1,num_it=1,mu=1.0,type="superbee",drag_x=0.0,drag_y=0.0,tau_margin=0.0,mask = false, beta=false)
    t=0
    Nx = grid.Nx
    if mask == false
        mask = ones(Nx+2)
    end
    for timeCount in 1:num_steps 
        #dt=update_vel_LW(vel,T,h,U_inflow,V_inflow,grid,dim,visc_prms,num_steps=num_it,mu=mu,type=type,drag_x=drag_x,drag_y=drag_y,tau_margin=tau_margin)
        #dt=update_vel(vel,T,h,U_inflow,V_inflow,grid,dim,num_steps=num_it,mu=mu,tau_margin=tau_margin)
        dt=update_vel(vel,T,h,U_inflow,V_inflow,grid,dim,visc_prms;num_steps=1,mu=1.0,tau_margin=tau_margin,mask=mask,beta=beta)
        update_thick(vel,T,h,U_inflow,V_inflow,H_inflow,accum,grid.dx,grid.dt*num_it,dim,type=type)
        # Gradient at the end point
        #h[:,end]=h[:,end-1]+(h[:,end-3]-4*h[:,end-2]+3*h[:,end-1])
        #T.xx[:,end]=T.xx[:,end-1]+(T.xx[:,end-3]-4*T.xx[:,end-2]+3*T.xx[:,end-1])
        #T.yy[:,end]=T.yy[:,end-1]+(T.yy[:,end-3]-4*T.yy[:,end-2]+3*T.yy[:,end-1])
        #h[:,end]= h[:,end-1]# .+ (h[:,end-1]-h[:,end-2])
        t+=grid.dt*num_it
    end
    return t
end

#=
function update_vel_flux_limiter(vel,T,h,U_inflow,V_inflow,grid,prms;,mu=1.0,type="superbee",drag_x=0.0,drag_y=0.0)
    u = vel.u 
    v = vel.v 

    txx = T.xx 
    tyy = T.yy 
    txy = T.xy

    
    dt = grid.dt 
    dx = grid.dx
    froude_number_squared = prms.froude_number_squared
    argand_upon_mach2 = prms.argand_upon_mach2
    mach_number_squared = prms.mach_number_squared
    froude_upon_mach2 = prms.froude_upon_mach2
    argand_number = prms.argand_number


    # Square and ratio of some non-dimensional parameters
    #froude_number_squared = froude_number^2
    #mach_number_squared = mach_number^2
    #froude_upon_mach2 = froude_number_squared/mach_number_squared
    #argand_upon_mach2 = argand_number / mach_number_squared


    thick_edgeFace = h
    thick = cellCenter(h)
    thick_cellCenter = thick
    #thick = ice_thick.thick_cellCenter
    #thick_grad_x,thick_grad_y = grad_thick(thick_cellCenter,grid)
    
    thick_sy =  0.5*(thick_cellCenter[:,2:end]+thick_cellCenter[:,1:end-1])
    thick_sx =  0.5*(thick_cellCenter[2:end,:]+thick_cellCenter[1:end-1,:])
    #thick_sx,thick_sy = cellCenter(ice_thick.thick_cellCenter)
    #h_end=0.5*(thick_cellCenter[1:end-1,end]+thick_cellCenter[2:end,end])

    set_bc(u,v,T.xx,T.yy,T.xy,thick,U_inflow,V_inflow)


    t=0
    σxx = T.xx + 0.5*T.yy
    σyy = T.yy + 0.5*T.xx
    σxy = 0.5*T.xy
    
    
    for timeCount in 1:num_steps
       
        η = mu.*visc(T,thick_edgeFace)

        # Calculate terms for later use
        σxx_minus_h2 = σxx- thick_edgeFace.^2
        σyy_minus_h2 = σyy- thick_edgeFace.^2

        # Calculate diffusion terms (hu_x)_x, (hu_y)_y
        hu_xx,hu_yy = diffusivity(thick_edgeFace,u,dx)
        # Calculate diffusion terms (hv_x)_x, (hv_y)_y
        hv_xx,hv_yy = diffusivity(thick_edgeFace,v,dx)

        hv_y = (v[3:end,:]-v[1:end-2,:]).*thick_edgeFace[2:end-1,:]/(2*dx)
        hu_y = (u[3:end,:]-u[1:end-2,:]).*thick_edgeFace[2:end-1,:]/(2*dx)
        hu_x = (u[:,3:end]-u[:,1:end-2]).*thick_edgeFace[:,2:end-1]/(2*dx)
        hv_x = (v[:,3:end]-v[:,1:end-2]).*thick_edgeFace[:,2:end-1]/(2*dx)

        
        σxy_y   = (σxy[3:end,:]-σxy[1:end-2,:])./thick_edgeFace[2:end-1,:]/(2*dx) 
        σxy_x   = (σxy[:,3:end]-σxy[:,1:end-2])./thick_edgeFace[:,2:end-1]/(2*dx)

        σxx_minus_h2_x = (σxx_minus_h2[:,3:end]-σxx_minus_h2[:,1:end-2])./thick_edgeFace[:,2:end-1]/(2*dx)
        σxx_minus_h2_x[:,1] .= (σxx_minus_h2[:,3]- σxx_minus_h2[:,2])/dx./thick_edgeFace[:,2]
        σxx_minus_h2_x[:,end] .= (σxx_minus_h2[:,end]- σxx_minus_h2[:,end-1])/dx./thick_edgeFace[:,end-1]


        σyy_minus_h2_y = (σyy_minus_h2[3:end,:]-σyy_minus_h2[1:end-2,:])./thick_edgeFace[2:end-1,:]/(2*dx)
        σyy_minus_h2_y[1,:] .= (σyy_minus_h2[3,:]- σyy_minus_h2[2,:])/dx./thick_edgeFace[2,:]# Assume gradient vanishes at inflow boundary . . .
     
        # Separate into waves traveling east/west and north/south
        u1 = u + 0.5*(froude_number_squared.*thick_edgeFace.*u/mach_number+σxx)
        u2 = u + 0.5*(-froude_number_squared.*thick_edgeFace.*u/mach_number+σxx)
        v1 = v + 0.5*(froude_number_squared.*thick_edgeFace.*v/mach_number+σyy)
        v2 = v + 0.5*(-froude_number_squared.*thick_edgeFace.*v/mach_number+σyy)

       

        source_1 =  -0.5*froude_number_squared.*u[2:end-1,2:end-1].*(thick_edgeFace[2:end-1,3:end].-thick_edgeFace[2:end-1,1:end-2])/dx/2/mach_number_squared
        source_2 =   0.5*froude_number_squared.*u[2:end-1,2:end-1].*(thick_edgeFace[2:end-1,3:end].-thick_edgeFace[2:end-1,1:end-2])/dx/2/mach_number_squared

        source_x[:,end] .= 0.0
        source_y[:,end] .= 0.0

        c1 = u .- 1/mach_number
        vel_seismic1 = Velocity(c1,c1)
        Fx1,junk = flux(vel_seismic1,w1,dx,dt;type=type)
        s1,junk = flux(vel_seismic1,0.25*froude_number_squared*thick_edgeFace.^2,dx,dt;type=type)
        c2 = u .+ 1/mach_number
        vel_seismic2 = Velocity(c2,c2)
        Fx2,junk = flux(vel_seismic2,w2,dx,dt;type=type)
        s2,junk = flux(vel_seismic2,-0.25*froude_number_squared*thick_edgeFace.^2,dx,dt;type=type)
        Fx1[:,end] .= Fx1[:,end-1]  
        s2[:,end] .= s2[:,end-1]
        s1[:,end] .= s1[:,end-1]

        w1[2:end-1,2:end-1] .= w1[2:end-1,2:end-1] .- dt/dx*(Fx1[:,2:end]-Fx1[:,1:end-1]) .+
                    source_x.*dt
                    #dt/dx.*(s1[:,2:end]-s1[:,1:end-1])

        w2[2:end-1,2:end-1] .= w2[2:end-1,2:end-1] .- dt/dx*(Fx2[:,2:end]-Fx2[:,1:end-1]) .+
                    source_y.*dt 
                    #dt/dx.*(s2[:,2:end]-s2[:,1:end-1])
           



        u[2:end-1,2:end-1] = u[2:end-1,2:end-1] +
            dt/froude_number_squared.*σxx_minus_h2_x[2:end-1,:] +
            # dt*σxy,y
            dt/froude_number_squared./thick_edgeFace[2:end-1,2:end-1].*(σxy[3:end,2:end-1]-σxy[1:end-2,2:end-1])/(2*dx) +
            # 0.5*dt^2*(hu_x)_x
            0.5*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*hu_xx[2:end-1,2:end-1] + 
            # 0.5*dt^2*0.5*(hv_y)_x  
            0.25*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hv_y[:,3:end]-hv_y[:,1:end-2])/(2*dx) +
            #-0.5*dt^2*Ar*(σxx/η)_x
            # Shear stress terms
            -0.5*dt^2*argand_upon_mach2.*(σxx[2:end-1,3:end]./η[2:end-1,3:end] - σxx[2:end-1,1:end-2]./η[2:end-1,1:end-2])/(2*dx)./thick_edgeFace[2:end-1,2:end-1] +
            # 0.5*dt^2*(h*u_y)_y 
            0.125*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*hu_yy[2:end-1,2:end-1] +   
            # 0.5*dt^2*(h*v_x)_y
            0.125*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hv_x[3:end,:]-hv_x[1:end-2,:])/(2*dx) +
            # 0.5*dt^2*2*Ar*(σxy/η)_y
            -0.25*(dt^2)*argand_upon_mach2./thick_edgeFace[2:end-1,2:end-1].*(σxy[3:end,2:end-1]./η[3:end,2:end-1]-σxy[1:end-2,2:end-1]./η[1:end-2,2:end-1])/(2*dx)
        u[1,:] = u[2,:]
        u[end,:]=u[end-1,:]
        u[:,end] = u[:,end-1]
        u[:,1]=2*U_inflow .- u[:,2] 
        
        
       
        v[2:end-1,2:end-1] = v[2:end-1,2:end-1] + 
            dt/froude_number_squared.*σyy_minus_h2_y[:,2:end-1] +
            dt/froude_number_squared./thick_edgeFace[2:end-1,2:end-1].*(σxy[2:end-1,3:end]-σxy[2:end-1,1:end-2])/(2*dx) +
            0.5*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*hv_yy[2:end-1,2:end-1] +        
            0.25*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hu_x[3:end,:]-hu_x[1:end-2,:])/(2*dx) -
            0.5*dt^2*argand_upon_mach2.*(σyy[3:end,2:end-1]./η[3:end,2:end-1] - σyy[1:end-2,2:end-1]./η[1:end-2,2:end-1])/(2*dx)./thick_edgeFace[2:end-1,2:end-1] +            
            # Shear stress terms
            0.125*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*hv_xx[2:end-1,2:end-1] +   
            0.125*dt^2/mach_number_squared./thick_edgeFace[2:end-1,2:end-1].*(hu_y[:,3:end]-hu_y[:,1:end-2])/(2*dx) +
            -0.25*dt^2*argand_upon_mach2./thick_edgeFace[2:end-1,2:end-1].*(σxy[2:end-1,3:end]./η[2:end-1,3:end]-σxy[2:end-1,1:end-2]./η[2:end-1,1:end-2])/(2*dx)
    
       

        v[1,:] = 2*V_inflow .- v[2,:]
        v[end,:]=2*V_inflow .- v[end-1,:]
        v[:,1]= v[:,2]       # Average across cells = margin vel
        v[:,end]= v[:,end-1] # Average across cells = margin vel
        
        σxx[2:end-1,2:end-1] = σxx[2:end-1,2:end-1] + 
            (dt*froude_upon_mach2 .- 0.5*dt^2*froude_upon_mach2^2*argand_number./η[2:end-1,2:end-1]).*(hu_x[2:end-1,:] .+0.5*hv_y[:,2:end-1]  - argand_number*σxx[2:end-1,2:end-1]./η[2:end-1,2:end-1]) +
            0.5*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxx_minus_h2[2:end-1,1:end-2]./thick_sx[:,1:end-1] .- (1 ./thick_sx[:,1:end-1] .+ 1 ./thick_sx[:,2:end]).*σxx_minus_h2[2:end-1,2:end-1] + σxx_minus_h2[2:end-1,3:end]./thick_sx[:,2:end])/(dx^2) +    
            0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σyy_minus_h2[1:end-2,2:end-1]./thick_sy[1:end-1,:] .- (1 ./thick_sy[1:end-1,:] .+ 1 ./thick_sy[2:end,:]).*σyy_minus_h2[2:end-1,2:end-1] + σyy_minus_h2[3:end,2:end-1]./thick_sy[2:end,:])/dx^2 +
            # Shear stress terms
            0.5*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy_y[:,3:end]-σxy_y[:,1:end-2])/(2*dx) +
            0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy_x[3:end,:]-σxy_x[1:end-2,:])/(2*dx)

        σxx[1,:]  = σxx[2,:]
        σxx[end,:] =σxx[end-1,:]
        σxx[2:end,end] = 2.0.*thick[:,end].^2  .- σxx[2:end,end-1]        
        σxx[:,1] = σxx[:,2]
    
        
        σyy[2:end-1,2:end-1] = σyy[2:end-1,2:end-1] + 
            (dt*froude_upon_mach2 .- 0.5*dt^2*froude_upon_mach2^2*argand_number./η[2:end-1,2:end-1]).*(hv_y[:,2:end-1] .+0.5*hu_x[2:end-1,:]  - argand_number*σyy[2:end-1,2:end-1]./η[2:end-1,2:end-1]) +
            0.5*dt^2/mach_number_squared*(σyy_minus_h2[1:end-2,2:end-1]./thick_sy[1:end-1,:] .- (1 ./thick_sy[1:end-1,:] .+ 1 ./thick_sy[2:end,:]).*σyy_minus_h2[2:end-1,2:end-1] + σyy_minus_h2[3:end,2:end-1]./thick_sy[2:end,:])/(dx^2).*thick_edgeFace[2:end-1,2:end-1] +    
            0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxx_minus_h2[2:end-1,1:end-2]./thick_sx[:,1:end-1] .- (1 ./thick_sx[:,1:end-1] .+ 1 ./thick_sx[:,2:end]).*σxx_minus_h2[2:end-1,2:end-1] + σxx_minus_h2[2:end-1,3:end]./thick_sx[:,2:end])/dx^2 +
            # Shear stress terms
            0.5*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy_x[3:end,:]-σxy_x[1:end-2,:])/(2*dx) +
            0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy_y[:,3:end]-σxy_y[:,1:end-2])/(2*dx)
        σyy[1,:] =  σyy[2,:]
        σyy[end,:] =  σyy[end-1,:]
        σyy[:,end]=  σyy[:,end-1]
        σyy[:,1]= σyy[:,2]
      
        
        
        # Double check?   
        σxy[2:end-1,2:end-1] = σxy[2:end-1,2:end-1] + 
            (0.5*dt*froude_upon_mach2 .- 0.25*dt^2*froude_upon_mach2^2*argand_number./η[2:end-1,2:end-1]).*(0.5*(hu_y[:,2:end-1] + hv_x[2:end-1,:]) - 2*argand_number*σxy[2:end-1,2:end-1]./η[2:end-1,2:end-1]) +
            0.125*dt^2/mach_number_squared.*(σxx_minus_h2_x[3:end,:]-σxx_minus_h2_x[1:end-2,:])/(2*dx) +
            0.125*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σyy_minus_h2_y[:,3:end]-σyy_minus_h2_y[:,1:end-2])/(2*dx) +
            0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy[2:end-1,1:end-2]./thick_sx[:,1:end-1] .- (1 ./thick_sx[:,1:end-1] .+ 1 ./thick_sx[:,2:end]).*σxy[2:end-1,2:end-1] + σxy[2:end-1,3:end]./thick_sx[:,2:end])/(dx^2) +
            0.25*dt^2/mach_number_squared*thick_edgeFace[2:end-1,2:end-1].*(σxy[1:end-2,2:end-1]./thick_sy[1:end-1,:] .- (1 ./thick_sy[1:end-1,:] .+ 1 ./thick_sy[2:end,:]).*σxy[2:end-1,2:end-1] + σxy[3:end,2:end-1]./thick_sy[2:end,:])/(dx^2)
        σxy[:,end] = -σxy[:,end-1];σxy[:,1]  = σxy[:,2]
        σxy[1,:]   = -σxy[2,:];    σxy[end,:]= -σxy[end-1,:] 
       
        σxy=σxy/2
       
        T.xx[1:end,1:end] =  (σxx-0.5*σyy)/0.75
        T.yy[1:end,1:end] =  (σyy-0.5*σxx)/0.75
        T.xy[1:end,1:end] = 2*σxy       
        t += dt  
             
    end
    return t
end
=#

