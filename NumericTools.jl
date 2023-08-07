#=
Module with some tools and utilities to do Lax Wendroff type integration
=#

#module NumericTools

#export advect, RectGrid
include("ice_constants.jl")
using Parameters
using LinearAlgebra


@with_kw mutable struct RectGrid
    Nx::Int64 = 100
    Lx::Float64 = 100e3
    Ly::Float64 = 25e3
    Ny::Int64 = floor(Int,Ly/Lx*Nx)
    dx::Float64 = 1/(Nx)
    xx::Vector{Float64} = LinRange(0,1,Nx+1)
    x::Vector{Float64} = LinRange(dx/2,1-dx/2,Nx)
    yy::Vector{Float64} = LinRange(0,Ly/Lx,Ny+1)
    y::Vector{Float64}= LinRange(dx/2,Ly/Lx-dx/2,Ny)
    dt::Float64 = 0.5
end

function flux_div(qx,qy,dx)
    qx_x   = (qx[2:end-1,3:end]-qx[2:end-1,1:end-2])/(2*dx)
    qy_y   = (qy[3:end,2:end-1]-qy[1:end-2,2:end-1])/(2*dx) 
    return qx_x+qy_y
end

function flux_div_staggered_grid(qx,qy,dx)
    qx_x   = (qx[:,2:end]-qx[:,1:end-1])/(dx)
    qy_y   = (qy[2:end,:]-qy[1:end-1,:])/(dx)
    qx_x =  0.5*(qx_x[1:end-1,:]+qx_x[2:end,:])
    qy_y =  0.5*(qy_y[:,1:end-1]+qy_y[:,2:end])
    return qx_x+qy_y
end


function advect(vel,h,accum,dx,dt; type="None")
    #=
    Advect quantity h using Lax-Wendroff
    Input: vel --velocity vel.u,vel.v
           h -- ice thickness or scalar variable to advect 
           accum  -- accumulation rate or right hand side 
           h_in -- ice thickness at grounding line or other BC
           grid -- grid structure
    =#

    #dt = grid.dt 
    #dx = grid.dx

    # Define fluxes
    qx = h.*vel.u
    qy = h.*vel.v
    

    # Flux divergence: (uh)_x + (vh)_y
    divHU = flux_div(qx,qy,dx)

    # Flux divergence of (accum*h)_x + (accum*h)_y
    divHa  = flux_div(accum.*vel.u,accum.*vel.v,dx)

    #accum_staggered_grid = cellCenter(accum)
    #thick_staggered_grid = cellCenter(h)

    # div(div(hu)u)
    #u_sx = 0.5*(vel.u[:,1:end-1] + vel.u[:,2:end])
    #u_staggered_grid = 0.5*(u_sx[1:end-1,:] + u_sx[2:end,:])
    u_staggered_grid = cellCenter(vel.u)
    #v_sy = 0.5*(vel.v[:,1:end-1] + vel.v[:,2:end])
    #v_staggered_grid = 0.5*(v_sy[1:end-1,:] + v_sy[2:end,:])
    v_staggered_grid = cellCenter(vel.v)
    divHU_staggered_grid = flux_div_staggered_grid(qx,qy,dx)
    #divHa = flux_div_staggered_grid(accum_staggered_grid.*u_staggered_grid,accum_staggered_grid.*v_staggered_grid,dx)
    qqx = divHU_staggered_grid.*u_staggered_grid
    qqy = divHU_staggered_grid.*v_staggered_grid
    div2HU = flux_div_staggered_grid(qqx,qqy,dx)
   

    # Update ice thickness based on Lax-Wendroff 2nd order update
    h[2:end-1,2:end-1] = h[2:end-1,2:end-1] .+ dt*(accum[2:end-1,2:end-1] .- divHU) .- 0.5*dt^2*(divHa .- div2HU)
    #h[:,1]= 2*h_in .- h[:,2]
    #h[:,end] .= h[:,end-1]
    #h[end,:] .= h[end-1,:]
    #h[1,:] .= h[2,:]
    #h[h.<0.1] .= 0.1
    return h
end

function pad_grid(u_cellCenter)

    Ny,Nx = size(u_cellCenter)
    Ny = Ny-1
    Nx = Nx-1

    
    # Initialize variables
    u = zeros(Ny+2,Nx+2)
    for i in 2:Ny+1
        for j in 2:Nx+1
            u[i,j]=0.25*(u_cellCenter[i-1,j-1]+u_cellCenter[i,j-1]+u_cellCenter[i-1,j]+u_cellCenter[i,j])
        end
    end
    u[:,1]=u[:,2]
    u[:,end]=u[:,end-1]
    u[1,:]=u[2,:]
    u[end,:]=u[end-1,:]
    return u
end

function cellCenter(var)
    # Calculate variable on staggered grid (sg)
    var_cellCenter =  0.5*(var[:,2:end]+var[:,1:end-1])
    var_cellCenter =  0.5*(var_cellCenter[2:end,:]+var_cellCenter[1:end-1,:])
    return var_cellCenter
end

function diffusivity(thick,var,dx)
    #=
    Calculate diffusivity terms
    diff_x = (D*u_x)_x
    diff_y = (D*u_y)_y
    =#

    # Preallocate an array to make sure it is the right size
    diff_x = zeros(size(thick))
    diff_y = zeros(size(thick))

    diff_x[2:end-1,2:end-1]=0.5*(var[2:end-1,3:end].*(thick[2:end-1,3:end]+thick[2:end-1,2:end-1]) 
        - var[2:end-1,2:end-1].*(thick[2:end-1,3:end]+2*thick[2:end-1,2:end-1]+thick[2:end-1,1:end-2])
        + var[2:end-1,1:end-2].*(thick[2:end-1,2:end-1]+thick[2:end-1,1:end-2]))/dx^2

    diff_y[2:end-1,2:end-1]=0.5*(var[3:end,2:end-1].*(thick[3:end,2:end-1]+thick[2:end-1,2:end-1]) +
        - var[2:end-1,2:end-1].*(thick[3:end,2:end-1]+2*thick[2:end-1,2:end-1]+thick[1:end-2,2:end-1])+
        + var[1:end-2,2:end-1].*(thick[2:end-1,2:end-1]+thick[1:end-2,2:end-1]))/dx^2
    
    return diff_x, diff_y
end

function advect_staggered_grid(vel,h,accum,dx,dt)
    accum_staggered_grid = cellCenter(accum)
    u_staggered_grid = cellCenter(vel.u)
    v_staggered_grid = cellCenter(vel.v)
    h_staggered_grid = cellCenter(h)
    vel_staggered_grid=Velocity(u_staggered_grid,v_staggered_grid)
    hnew=advect(vel_staggered_grid,h_staggered_grid,accum_staggered_grid,dx,dt)
    #h[2:end-1,2:end-1]=cellCenter(hnew)
    return hnew
end


function limiter(r ; type="superbee")
    if type == "vanlear"
        phi = (r .+ abs.(r))./(1 .+ abs.(r))
    elseif type == "superbee"
        phi = max.(0, min.(1, 2*r), min.(2, r))
    elseif type == "MC"
        phi = max.(0, min.((1 + r)./2, 2, 2*r))
    elseif type =="Laxwendroff"
        phi =1.0
    else
        phi = 0.0
    end

    return phi
end

function flux(vel,q,dx,dt;type="superbee")
    
    #dt = grid.dt 
    #dx = grid.dx
    # Calculate velocity at mid-point of cells
    umid_x=0.5*(vel.u[2:end-1,1:end-1]+vel.u[2:end-1,2:end])
    #vmid_x=0.5*(vel.v[2:end-1,1:end-1]+vel.v[2:end-1,2:end])

    #umid_y=0.5*(vel.u[1:end-1,2:end-1]+vel.u[2:end,2:end-1])
    vmid_y=0.5*(vel.v[1:end-1,2:end-1]+vel.v[2:end,2:end-1])

    n,m=size(umid_x)
    r1 = [ones(n,1) (q[2:end-1,2:end-1]-q[2:end-1,1:end-2])./(q[2:end-1,3:end]-q[2:end-1,2:end-1].+1e-12)]
    r2 = [(q[2:end-1,3:end]-q[2:end-1,2:end-1])./(q[2:end-1,2:end-1]-q[2:end-1,1:end-2].+1e-12) ones(n,1)]
    rx = r1.*(umid_x.>0) .+ r2.*(umid_x.<=0.0)
    #if any(isnan,rx)
    #    print("NAN's rx!!!")
    #end
    n,m=size(vmid_y)
    r1 = [ones(1,m) ; (q[2:end-1,2:end-1]-q[1:end-2,2:end-1])./(q[3:end,2:end-1]-q[2:end-1,2:end-1].+1e-12)]
    r2 = [(q[3:end,2:end-1]-q[2:end-1,2:end-1])./(q[2:end-1,2:end-1]-q[1:end-2,2:end-1].+1e-12) ; ones(1,m)]
    ry = r1.*(vmid_y.>0) .+ r2.*(vmid_y.<=0.0)
    #if any(isnan,rx)
    #    print("NAN's ry!!!")
    #end

    #r1 = (q[2:end-1,2:end-1]-q[2:end-1,1:end-2])./(q[2:end-1,3:end]-q[2:end-1,2:end-1])
    #r2 = (q[2:end-1,3:end]-q[2:end-1,2:end-1])./(q[2:end-1,2:end-1]-q[3:end,2:end-1])
    #ry = [ones(n,1) r1.*(umid_x[:,1:end-1].>0) .+ r2.*(umid_x[:,2:end].<=0)]

   
    phix = limiter(rx,type=type)
    phiy = limiter(ry,type=type)
    #print(size(phi))
    #print(size(umid_x))
    # Flux at i-1/2
    #qq = q[2:end-1,1:end-1].*(umid_x.>0) .+ q[2:end-1,2:end].*(umid_x.<=0)
    Fx = 0.5*umid_x.*(q[2:end-1,1:end-1].*(umid_x.>0) .+ q[2:end-1,2:end].*(umid_x.<=0)) .+
         0.25*abs.(umid_x).*(1.0 .- abs.(umid_x*dt/dx)).*phix.*(q[2:end-1,2:end]-q[2:end-1,1:end-1])

    Fy = 0.5*vmid_y.*(q[1:end-1,2:end-1].*(vmid_y.>0) .+ q[2:end,2:end-1].*(vmid_y.<=0)) .+
         0.25*abs.(vmid_y).*(1.0 .- abs.(vmid_y*dt/dx)).*phiy.*(q[2:end,2:end-1]-q[1:end-1,2:end-1])
        
   return Fx,Fy
end

function advect_LW(vel,h,accum,dx,dt;type="superbee")
    Fx,Fy = flux(vel,h,dx,dt;type=type)
    #Fx[:,end] = Fx[:,end-1]
    h[2:end-1,2:end-1] = h[2:end-1,2:end-1] .+ 
            -dt/dx*(Fx[:,2:end]-Fx[:,1:end-1]) .+
            -dt/dx*(Fy[2:end,:]-Fy[1:end-1,:]) .+
            +accum[2:end-1,2:end-1].*dt
    return h
end

function principal_stess(grid,T,h,k)
    x = grid.x[1:k:end]
    y = grid.y[1:k:end]
    xx = zeros(length(x)*length(y))
    yy = zeros(length(x)*length(y))
    u1 = zeros(length(x)*length(y))
    v1 = zeros(length(x)*length(y))
    u2 = zeros(length(x)*length(y))
    v2 = zeros(length(x)*length(y))
    e1 = zeros(length(x)*length(y))
    e2 = zeros(length(x)*length(y))
    q = 1
    for i in k:k:length(grid.x)-k
        for j in k:k:length(grid.y)-k
            A = [T.xx[j,i]./h[j,i] T.xy[j,i]./h[j,i];T.xy[j,i]./h[j,i] T.yy[j,i]./h[j,i]]
            (evals, evecs) = eigen(Symmetric(A))
            u1[q] = evals[1]*evecs[1,1]
            u2[q] = evals[2]*evecs[2,1]
            v1[q] = evals[1]*evecs[1,2]
            v2[q] = evals[2]*evecs[2,2]
            xx[q] = grid.x[i]
            yy[q] = grid.y[j]
            e1[q] = evals[1]
            e2[q] = evals[2]
            q = q+1
        end
    end
    return xx,yy,u1,v1,u2,v2,e1,e2
end
