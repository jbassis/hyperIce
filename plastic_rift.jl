#=
Simulate plastic rifting and deformation using m-ice rheology
=#

# To do
# Save file as something ???
# Add code to specify margin plastic stress
# Add code to simulate evolution fo plastic strain (advect log of values?)
# See if we can simulate a rift??
# Step #1 Initialize ice shelf to steady-state
# Step #2 Allow plastic strain to accumulate

include("HyperIceShelf2D.jl")
using Plots; pyplot()
using Printf
import .ShelfModel2D as ice
using JLD2


#--------------------
# Define material properities and dimensionless numbers
# Material properties and stuff
mat = ice.Material()
mat.B = 3.2e8
c = sqrt(4*1e9/mat.ρ_i)
H0 = 1e3
Lx = L0 = 80e3  # length of ice shelf
accum =0.0/mat.seconds_in_year
H0 = 1.4e3 # Grounding line ice thickness
U_inflow = 1000/mat.seconds_in_year;   # Inflow velocity
U_peak = U_inflow*mat.seconds_in_year
V_inflow = 0.0                        # Margin no-penetration bc
#U0 = U_inflow
U0 = U00=(mat.di*mat.ρ_i*mat.g*H0/(4*mat.B))^mat.n*Lx;  # Ice velocity scale [m/s]
argand_number = (mat.di*mat.ρ_i*mat.g*H0/(4*mat.B))^mat.n*(Lx/U0)
T0 = (mat.di*mat.ρ_i*mat.g*H0/(4)) # Dimensionless stress



# Define dimensionless numbers
mach_number = 0.1
froude_number = 1.0
argand_number = 1.0
#froude_number =1.0#3.75
#mach_number = 0.1#0.05
U00 = U0 = (mat.di*mat.ρ_i*mat.g*H0/(4*mat.B))^mat.n*(Lx/argand_number)

#U00 = 30.0

#=
# Testing elastic limit
U00= 30.0
c = sqrt(4*1e9/mat.ρ_i)
argand_number = (mat.di*mat.ρ_i*mat.g*H0/(4*mat.B))^mat.n*(Lx/U00)
mach_number = U00/c 
# This ensures that the ratio of the Froude to Mach number remains appropriate
froude_upon_mach = c*sqrt(2)/sqrt(g*di*H0)
froude_number = froude_upon_mach*mach_number
=#

# End testing
dim=ice.Dimensionless(mach_number=mach_number,froude_number=froude_number,argand_number=argand_number)

#--------------------
# Create grid
# Numerical and geometric parameters
Ly = Lx/2 # width of ice shelf
Nx = 200 # number or grid points
grid = ice.RectGrid(Nx=Nx,Lx=Lx,Ly=Ly)
dx=grid.dx
dt = mach_number*dx
grid.dt = dt/3
#--------------------
# Define analytic initial condition for confined ice tongue
ha,ua,Exx=ice.steady_state(grid.xx,grid.Lx,accum,H0,U_inflow,mat)

# Initialize veriables based on analytic solution
vel, T, h, h_cellCenter=ice.initialize2D(grid,ua/U0,ha/H0,H0/H0)
ice.set_bc(vel.u,vel.v,T.xx,T.yy,T.xy,h,U_inflow/U0,V_inflow/U0)
uaa = 0.5*(ua[1:end-1]+ua[2:end])
haa = 0.5*(ha[1:end-1]+ha[2:end])
txx_old = copy(T.xx)
tyy_old = copy(T.yy)

# Make accumulation an array with the right size
accum =0.0/mat.seconds_in_year*ones(size(h))
#accum[1:3,:].=-1.0/mat.seconds_in_year
#accum[end-3:end,:].=-1.0/mat.seconds_in_year

#accum[31:71,:] .= -1.0/mat.seconds_in_year



# Create plastic strain array
epsII = ones(size(h))*1e-16

# Initialize viscosity parameters
visc_prms = ice.visc_prms(eps=epsII) 
tau_y =0.75#0.78#0.804
visc_prms.tau_y = 250e3/T0#tau_y*0.95#-0.03#*0.95#*0.83
visc_prms.tau_min = 0.1/2#0.04
visc_prms.eps_crit =0.1/8
visc_prms.alpha = 0.0
visc_prms.μ = 0.4
# Initialize rift
#start = Int(Nx*3/4)
#last = start + Int(round(3e3/(grid.dx*Lx)))
#eps[1:end,start:last].= 0.01
t = 0
num_it = 1
num_steps=1
global h_old = copy(h)



y = [0;grid.y*Lx/Ly;1]


slip_frac = 0.7
tau_margin = 175e3/T0
beta0 = 0.0


#fname = "../data/rift_initial_$Nx"*"_U_inflow_$U_peak"*"slip_frac_$slip_frac"*"yield_$tau_y"*".jld"
#fname = "../data/rift_initial_margin$Nx"*"_U_inflow_$U_peak"*"slip_frac_$slip_frac"*"yield_$tau_y"*".jld"
#fname = "../data/rift_initial_$Nx"*"_U_inflow_500.0"*"slip_frac_$slip_frac"*"yield_$tau_y"*".jld"
#fname = "../data/rift_initial_$Nx"*"_U_inflow_360.0"*"slip_frac_$slip_frac"*"yield_$tau_y"*".jld"
fname = "../data/plume/plume_initial_$Nx"*"_U_inflow_$U_peak"*"slip_frac_$slip_frac"*"yield_$tau_y"*".jld"
#fname = "../data/plume/rift_margin$Nx"*"_U_inflow_$U_peak"*"slip_frac_$slip_frac"*"yield_$tau_y"*".jld"
#fname = "../data/initial/rift_initial_$Nx"*"_U_inflow_$U_peak"*"slip_frac_$slip_frac"*"yield_$tau_y"*".jld"
fname = "../data/initial/rift_initial_$Nx"*"_U_inflow_$U_peak"*"slip_frac_$slip_frac"*"yield_$tau_y"*"margin_$tau_margin"*".jld"
#fname = "../data/initial/rift_initial_$Nx"*"_Lx_$Lx"*"_U_inflow_$U_peak"*"slip_frac_$slip_frac"*"yield_$tau_y"*"margin_$tau_margin"*".jld"
#fname = "../data/initial/rift_initial_friction_$beta0"*"_Nx_$Nx"*"_Lx_$Lx"*"_U_inflow_$U_peak"*"slip_frac_$slip_frac"*"yield_$tau_y"*"margin_$tau_margin"*".jld"
#tau_margin = 250e3/T0
#fname = "../data/margin/rift_initial_friction_0.0_Nx_200_Lx_80000.0_U_inflow_1000.0slip_frac_0.7yield_0.75margin_0.5027335821269093.jld"
d=JLD2.load(fname)
T.xx = d["Txx"]
T.yy = d["Tyy"]
T.xy = d["Txy"]
h = d["h"];h_old = copy(h)
vel.u = d["u"]*U0/U00
vel.v = d["v"]*U0/U00
U0=U00

#accum = -d["m"]
# 200 kPa
#T0 = (mat.di*mat.ρ_i*mat.g*H0/(4))
#tau_y = 250e3/T0
visc_prms.tau_y = 40e3#250e3/T0
#tau_margin = 200e3/T0

E0 = zeros(size(T.xx))
#slip_frac = 0.7
pow = 2.0/3
U_in = 2^(pow-1)*2^(pow+1)*U_inflow*(1-slip_frac).*((y.*(1.0 .- y))).^pow .+ slip_frac*U_inflow
U_in = U_in
plastic_fail = true
De = mach_number^2/argand_number/froude_number^2
println("Debrah Number ",De," Dissipation Number ",mach_number^2/De^2," time step size ",dt*L0/U00)
mask = ones(Nx+2)
beta = zeros(grid.Ny,grid.Nx)
#beta0 = 0.0
#beta[Integer(grid.Ny/2)-1:Integer(grid.Ny/2)+1,133-1:133+1].= beta0
f=ones(Nx+2)#2*rand(Nx+2)
p=rand(Nx+2)
#mask[end-20:end].= 1.0
#prob = 0.1
#f[p.<prob].=0.0
#prob = 0.05
#f[p.<prob].=2.0

#f[end-50:end-48].=2.0


#mask[end-20:end].=0.
#mask[40:end].=0.

t = 0.0
for i in 0:1000
    tt = t*Lx/U0/mat.seconds_in_year
    #if tt>200
    #    mask[40:end].= max(1.0 - (tt-200.0)/10,0.0)
    #end
    #global g = min(t*Lx/U0/mat.seconds_in_year/2,4)
    #global U_in =  U_inflow*(g+1)
    #println("Maximum of uin ",maximum(U_in*mat.seconds_in_year)," ",t*Lx/U0/mat.seconds_in_year/2," g ",g+1)
    #if plastic_fail==true 
    #    if mod(i,2)==0
    #        local fname = "../data/rift/ex_tabular/rift_margin$Nx"*"_U_inflow_$U_peak"*"slip_frac_$slip_frac"*"yield_$tau_y"*"_"*string(i)*".jld"
    #        println(fname)
    #        JLD2.save(fname,"Txx", T.xx, "Tyy", T.yy, "Txy", T.xy, "u", vel.u, "v", vel.v, "h", h,"x",grid.x,"y",grid.y,"t",t*Lx/U0/mat.seconds_in_year,"epsII",epsII)
    #    end
    #end
    for j in 1:200
        # Update 
        local dt=ice.update(vel,T,h,U_in/U0,V_inflow/U0,H0/H0,accum*Lx/U0/H0,grid,dim,visc_prms;num_steps=num_steps,num_it=num_it,type="superbee",tau_margin=tau_margin.*mask.*f,mask=mask,beta=beta)
        #local dt=ice.update_vel(vel,T,h,U_in/U0,V_inflow/U0,grid,dim,visc_prms,num_steps=1,tau_margin=tau_margin) 
        
        if plastic_fail == true
            global epsII= ice.update_strain(T,h,visc_prms,dt)
            #epsII = min.(epsII,visc_prms.eps_crit)
            epsII = max.(epsII,1e-16)
            local q = log10.(epsII).*h       
            ice.advect_LW(vel,q,q.*accum*Lx/U0/H0,grid.dx,grid.dt,type="superbee")
            epsII[:,:] .= 10.0.^(q./h)
            epsII[:,1] .=  .- epsII[:,2]
            epsII[:,end] .= epsII[:,end-1]
            epsII[end,:] .= epsII[end-1,:]
            epsII[1,:] .= epsII[2,:]
            epsII = min.(epsII,visc_prms.eps_crit)
            epsII = max.(epsII,1e-16)
            epsII[:,1:10] .= 1e-16
            visc_prms.eps = copy(epsII)
        end
        
        global t += dt
    end
    if any(isnan,h)==true 
        break 
    end

    if any(isnan,vel.u)==true 
        break 
    end
    
    global err = sum(sqrt.((h_old-h).^2))/(grid.Nx*grid.Ny)*H0
    #local err = sum(sqrt.((v_old-vel.u).^2))/(grid.Nx*grid.Ny)*U0*mat.seconds_in_year
    global h_old = copy(h)
    println(i," time ",t*Lx/U0/mat.seconds_in_year, " L2 norm ",err)
    #if err<1e-3
    #    break
    #end
    #if plastic_fail==false
        #local fname = "../data/rift_initial_$Nx"*"slip_frac_$slip_frac"*"yield_$tau_y"*".jld"
        #U_peak = U_inflow*mat.seconds_in_year
        #local fname = "../data/rift_initial_margin$Nx"*"_U_inflow_$U_peak"*"slip_frac_$slip_frac"*"yield_$tau_y"*".jld"
        #local fname = "../data/plume/rift_margin$Nx"*"_U_inflow_$U_peak"*"slip_frac_$slip_frac"*"yield_$tau_y"*".jld"
        #JLD2.save(fname,"Txx", T.xx, "Tyy", T.yy, "Txy", T.xy, "u", vel.u, "v", vel.v, "h", h,"x",grid.x,"y",grid.y)
    #end
    

    #=
    local plot1=plot(grid.x*Lx/1e3,vel.u[floor(Int,grid.Ny/2),2:end-1]*U0,xlabel="distance [km]",ylabel="velocity [dimensionless]",legend=false,linewidth=2);
    plot!(grid.xx*Lx/1e3,ua,color=:red,linestyle=:dash,xlabel="distance [km]",ylabel="velocity [dimensionless]",legend=false);
    #plot!(ylims=(0.2,0.6))
    plot!(xlims=(0,Lx/1e3))
    #local plot3=plot(grid.x*Lx/1e3,mu[floor(Int,grid.Ny/2),2:end-1],xlabel="distance [km]",ylabel="mu",legend=false,linewidth=2);
    #plot!(grid.x*Lx/1e3,mu_old[floor(Int,grid.Ny/2),2:end-1],xlabel="distance [km]",ylabel="mu",legend=false,linewidth=2,color=:red);

    #plot!(ylims=(0,1/visc_ratio))
    plot2=plot(grid.x*Lx/1e3,h[floor(Int,grid.Ny/2),2:end-1],xlabel="distance [km]",ylabel="thickness [dimensionless]",legend=false,linewidth=2);
    plot!(grid.x*Lx/1e3,haa/H0,color=:red,linestyle=:dash,xlabel="distance [km]",ylabel="thickness [dimensionless]",legend=false,linewidth=2);
    plot!(xlims=(0,Lx/1e3))
    #plot!(ylims=(0,1.0./visc_ratio))
    
    plot(plot1, plot2, layout=(2,1))
    
    #plot1=contourf(vel.u[2:end-1,2:end-1])
    =#
    a = (T.xx .+ T.yy)
    b = (T.xx .- T.yy)
    global T1 = a .+ 0.5*sqrt.(b.^2 +T.xy.^2)
    global T2 = a .- 0.5*sqrt.(b.^2 +T.xy.^2)
    
    global E1 = (T1.^2).*(T1.>0) .+ (T2.^2).*(T2.>0)
    global dE = E1-E0
    global speed= sqrt.(vel.u.^2+vel.v.^2)
    local plot2=contourf(grid.x,grid.y,log10.(epsII[2:end-1,2:end-1]),aspect_ratio=:equal)
    
    local plot1=contourf(grid.x,grid.y,h[2:end-1,2:end-1],c=:viridis,aspect_ratio=:equal,clims=(0.55,1.0))
    global tau_e = sqrt.(T.xx.^2 + T.yy.^2 + T.xy.^2 + T.xx.*T.yy)
    s = 0.1
    local xx,yy,u1,v1,u2,v2,e1,e2 = ice.principal_stess(grid,T,h,20)

    quiver!(xx,yy,quiver=(-s*u2.*(e2.>0),-s*v2.*(e2.>0)),color=:black,linewidth=:1,aspect_ratio=:equal)
    quiver!(xx,yy,quiver=(s*u2.*(e2.>0),s*v2.*(e2.>0)),color=:black,linewidth=:1,aspect_ratio=:equal)
    quiver!(xx,yy,quiver=(-s*u1.*(e1.>0),-s*v1.*(e1.>0)),color=:black,linewidth=:1,aspect_ratio=:equal)
    quiver!(xx,yy,quiver=(s*u1.*(e1.>0),s*v1.*(e1.>0)),color=:black,linewidth=:1,aspect_ratio=:equal)
    
    quiver!(xx,yy,quiver=(-s*u2.*(e2.<0),-s*v2.*(e2.<0)),color=:red,linewidth=:1,aspect_ratio=:equal)
    quiver!(xx,yy,quiver=(s*u2.*(e2.<0),s*v2.*(e2.<0)),color=:red,linewidth=:1,aspect_ratio=:equal)
    quiver!(xx,yy,quiver=(-s*u1.*(e1.<0),-s*v1.*(e1.<0)),color=:red,linewidth=:1,aspect_ratio=:equal)
    quiver!(xx,yy,quiver=(s*u1.*(e1.<0),s*v1.*(e1.<0)),color=:red,linewidth=:1,aspect_ratio=:equal)
   
    local plot3=contourf(grid.x,grid.y,speed[2:end-1,2:end-1]*U00*mat.seconds_in_year,c=:viridis,aspect_ratio=:equal)

    #filter1 = e1>0
    #filter2 = e2>0

    
    #local plot3=contourf(grid.x,grid.y,dE[2:end-1,2:end-1],c=:viridis,aspect_ratio=:equal)
    #local plot2=contourf(grid.x,grid.y,tau_e[2:end-1,2:end-1]*T0/1e3,c=:viridis,aspect_ratio=:equal,clims=(100,350))
    #plot(plot3,plot1, plot2, layout=(3,1))
    plot(plot3,plot1,plot2, layout=(3,1))

    #display()
    gui()
    
    if err<1e-3
        break
    end
    global E0=E1

end
time=@sprintf("time %2.2f a",t*Lx/U0/mat.seconds_in_year)
println(time)
       

if plastic_fail==false
    U_peak = U_inflow*mat.seconds_in_year
    fname = "../data/margin/rift_initial_friction_$beta0"*"_Nx_$Nx"*"_Lx_$Lx"*"_U_inflow_$U_peak"*"slip_frac_$slip_frac"*"yield_$tau_y"*"margin_$tau_margin"*".jld"
    JLD2.save(fname,"Txx", T.xx, "Tyy", T.yy, "Txy", T.xy, "u", vel.u, "v", vel.v, "h", h,"x",grid.x,"y",grid.y,"beta",beta,"mask",mask)
end





