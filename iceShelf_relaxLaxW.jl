using Plots; pyplot()
using Printf
using JLD2
using LaTeXStrings

include("SteadyStates.jl")

# Define which experiment
experiment = 3

# Physical constants
B0 = 3.2e8
seconds_in_year = 365.2422*24*60*60
n = 3.0;                            # Flow-law exponent for ice [non-dimensional]% Conversion factors
rho_ice = 910;                      # Density of ice (assumed constant) [kg/m^3]
rho_w = 1028;                       # Density of sea water (assumed constant) [kg/m^3]
di = (1-rho_ice/rho_w)
g = 9.81;                           # Acceleration of gravity [m/s^2]
G = 1e9 # Shear modulus of elasticity

# Numerical and geometric parameters
L0 = 80e3
Nx = 320
dx = 1/(Nx)
xx = LinRange(0,1,Nx+1)
x = LinRange(dx/2,1-dx/2,Nx);

# Ice shelf characteristics and scales
accum = 0.0/seconds_in_year
H0 = 1.4e3
U0 = U00= (di*rho_ice*g*H0/(4*B0))^n*L0;  # Ice velocity scale [m/s]
U_in = U_inflow = 1000/seconds_in_year/U0;
P0 = B0*(U0/L0)^(1/n)/1e3;


# Initial condition
h,u,Exx=steady_state_old(xx,L0,0,B0,H0,U_inflow*U0)

# Scale terms
thick = h/H0
u     = u/U0
ha = thick # Analytic
thick_staggered_grid = 0.5*(thick[1:end-1] + thick[2:end])
ua = 0.5*(u[1:end-1] + u[2:end])
Exx   = Exx/(U0/L0)
p = thick.^2
pa = 0.5*(p[1:end-1] + p[2:end]);
Q0 = u[1]*thick[1];


#__________________________________________________________________________
# Example 0:  Test case
# Set non-dimensional parameters
mach_number    = 0.1;  
froude_number  = sqrt(10);
argand_number  = 1.0;

# Initial condition for elocity and ice thickness
velocity = [ua[1];ua;ua[end]];
velocity_new = copy(velocity);
stress =   [pa[1];pa;pa[end]];
stress_new = copy(stress);
U = 0

#__________________________________________________________________________

#__________________________________________________________________________
# Example 0:  Test that numerical solution matches analytic steady state solution
# Set non-dimensional parameters
if experiment==0
    mach_number    = 0.1;  
    froude_number  = 1.0;
    argand_number  = 1.0;
    up = copy(ua)
end

#__________________________________________________________________________
# Example 1  Numerical example 1 testing analytic perturbation
if experiment==1
    mach_number    = 1.0  
    froude_number  = 1.0
    argand_number  = 1.0
    
    print()
    #Velocity perturbation
    up=ua.*(1 .+ 0.2*exp.(-(x .- x[end]).^2/0.05^2))
end

if experiment==2
    mach_number    = 0.1  
    froude_number  = 1.0
    argand_number  = 1.0
    #Velocity perturbation
    up=ua.*(1 .+ 0.2*exp.(-(x .- x[end]).^2/0.05^2))
end


#__________________________________________________________________________
#Example 3:  Elastic
# Set non-dimensional parameters
if experiment==3
    # Define elastic wave speed and some physical quantities
    c=sqrt(4*G/rho_ice)
    froude_upon_mach = c*sqrt(2)/sqrt(g*di*H0)
    # Define velocity scale
    U00 = 30.0
    
    # Compute non-dimensional numbers
    argand_number = (di*rho_ice*g*H0/(4*B0))^n*(L0/U00)
    mach_number = U00/c
    froude_number = froude_upon_mach*mach_number

    
    ua = ua*U0/U00
    up = copy(ua)
    U_in = U_inflow*U0/U00
    De = mach_number^2/(argand_number*froude_number^2)
    println("Debrah Number ",De," Dissipation Number ",De^2/mach_number^2)

    # Scale perturbation velocity to relevant units
    U = 25e3/seconds_in_year/U00 # 25 km/a/U00

    #data=load("../Data/elastic_steady_state.jld")
    #ua = data["velocity"][2:end-1]
    #up = data["velocity"][2:end-1]
    #pa = data["stress"][2:end-1]
end

#Initial condition velocity and ice thickness
velocity = [up[1];up;up[end]];velocity_new = velocity;
stress =   [pa[1];pa;pa[end]];stress_new = stress;

# Determine some quantities that we will need
thick_edgeFace = 0.5*(thick[1:end-1] + thick[2:end]);
grad_thick = [0;(thick[3:end] - thick[1:end-2]);thick[end-2]-4*thick[end-1]+3*thick[end]]/2/dx;
#grad_thick = [0;(thick[3:end] - thick[1:end-2]);0.0]/2/dx;



# Square and ratio of some non-dimensional parameters
froude_number_squared = froude_number^2
mach_number_squared = mach_number^2
froude_upon_mach2 = froude_number_squared/mach_number_squared
argand_upon_mach2 = argand_number / mach_number_squared


timeCount=0;

time=@sprintf("Number of time steps: %2.2f",timeCount/Nx);




dt = mach_number*dx/2;
t=0
#__________________________________________________________________________


# Begin time-stepping
timeGap = 1
num_time_steps = Nx*100
seismogram = zeros(num_time_steps,4)
velocity_array = zeros(Nx,Integer(num_time_steps/Nx))
stress_array = zeros(Nx,Integer(num_time_steps/Nx))

#pulseWidth = 1200*3

t=0
pos = zeros(size(pa));

pulseWidth = 93556.34985743872*U00/L0/2#1200#30*60*U00/L0
pulseWidth = 90*U00/L0#*time_factor#*L0/80e3
#pulseWidth = 1200

De = mach_number^2/argand_number/froude_number_squared
println("Debrah Number ",De," Dissipation Number ",mach_number^2/De^2," time step size ",dt*L0/U00)
for timeCount in 1:num_time_steps
   

    step = (timeCount)/(2*Nx)
    if mod(timeCount-1,Nx)==0 && experiment<3
        #println(timeCount-1," ",(timeCount-1)/Nx)
        global step_count = Integer((timeCount-1)/(Nx))+1
        global velocity_array[:,step_count]=(velocity[2:end-1]-ua)./ua
        global stress_array[:,step_count]=(stress[2:end-1]-pa)./pa
    end


       
    
    if mod(timeCount,100) == 0
        #println(De/maximum(viscosity)," ",De," ", De/minimum(viscosity))
        #println(step)
        local time=@sprintf("Number of time steps: %2.2f",timeCount/Nx);
        #println(time)
        global plot1=plot(x*L0/1e3,((velocity[2:end-1])-ua)./ua,color=:black,linewidth=2,
            xlabel="distance [km]",ylabel="velocity [dimensionless]",
            legend=false,title=time)
        #global plot3=plot!(x*L0/1e3,ua,color=:red,linestyle=:dash,linewidth=2,
        #    xlabel="distance [km]",ylabel="analytic",
        #    legend=false,title=time)
        plot!(xlims=(0,L0/1e3))
        #plot!(ylims=(-1.5,1.5))
        #scatter!([0.0],[U_inflow],color = "red",markersize = 5)
        global plot2=plot(x*L0/1e3,((stress[2:end-1])-pa)./pa,color=:black,linewidth=2,
            xlabel="distance [km]",ylabel="stress [dimensionless]",
            legend=false)
        #global plot4=plot!(x*L0/1e3,pa,color=:red,linestyle=:dash,linewidth=2,
        #    xlabel="distance [km]",ylabel="analytic",
        #    legend=false,title=time)
        plot!(xlims=(0,L0/1e3))
        #plot!(ylims=(-1.5,1.5))
        #plot!(ylims=(-1,1))
        plot(plot1, plot2, layout=(1,2))
        gui()
    end
    
    
    
    
    global U_inflow = U_in + U*exp(-((t-10*pulseWidth)./(pulseWidth)).^2)
    global t = t+dt
    global velocity[[1,end]] = [(2*U_inflow-velocity[2]), velocity[end-1]]
    global stress[[1,end]] = [stress[2],2*thick[end]^2-stress[end-1]]
   
    # Define viscosity
    global viscosity = ((abs.(stress[2:end-1])./thick_edgeFace) .+ 1e-12).^(1-n);
    global viscosity = [viscosity[1];viscosity;viscosity[end]];
    
    global stress_minus_h2 = [stress[1]-thick_edgeFace[1].^2; stress[2:end-1] .- thick_edgeFace.^2; stress[end]-thick_edgeFace[end].^2]
    global stress_minus_h2[[1,end]] = [stress_minus_h2[2],stress_minus_h2[end-1]]

    global velocity[2:end-1] = velocity[2:end-1] + 
         0.5*(dt/dx)./(thick_edgeFace)./(froude_number_squared).*(stress_minus_h2[3:end]-stress_minus_h2[1:end-2]) +
         -0.25*(dt^2/dx)./(thick_edgeFace) * argand_upon_mach2.*(stress[3:end]./viscosity[3:end]-stress[1:end-2]./viscosity[1:end-2]) + 
         0.5*(dt/dx)^2/mach_number_squared.*(thick[1:end-1].*velocity[1:end-2] - (thick[1:end-1] + thick[2:end]).*velocity[2:end-1] + thick[2:end].*velocity[3:end])./thick_edgeFace#  +
         #-2*dt/froude_number_squared*0.5.*(grad_thick[1:end-1] + grad_thick[2:end]); 
    global velocity[[1,end]] = [(2*U_inflow-velocity[2]), velocity[end-1]]

    global stress[2:end-1] = stress[2:end-1] + 
        +dt*froude_upon_mach2*(0.5*thick_edgeFace.*(velocity[3:end]-velocity[1:end-2])/dx  - argand_number*stress[2:end-1]./viscosity[2:end-1]) +
        +0.5*(dt/dx)^2/mach_number_squared*(stress[1:end-2]./thick[1:end-1] .- (1 ./thick[1:end-1] .+ 1 ./thick[2:end]).*stress[2:end-1] + stress[3:end]./thick[2:end]).*thick_edgeFace +
        -dt^2/mach_number_squared*thick_edgeFace.*(grad_thick[2:end]-grad_thick[1:end-1])/dx +
        -0.5*dt^2*froude_upon_mach2^2*argand_number./viscosity[2:end-1].*(0.5*thick_edgeFace.*(velocity[3:end]-velocity[1:end-2])/dx  - argand_number*stress[2:end-1]./viscosity[2:end-1])
    global stress[[1,end]] = [stress[2],2*thick[end]^2-stress[end-1]]
    #global stress[[1,end]] = [stress[2],stress[end-1]]

    #=
    global velocity_new[2:end-1] = velocity[2:end-1] +
        + 0.5*(dt/dx)./(thick_edgeFace)./(froude_number_squared).*(stress[3:end]-stress[1:end-2]) +
        - 0.25*(dt^2/dx)./(thick_edgeFace) * argand_upon_mach2.*(stress[3:end]./viscosity[3:end]-stress[1:end-2]./viscosity[1:end-2]) +
        + 0.5*(dt/dx)^2/mach_number_squared.*(thick[1:end-1].*velocity[1:end-2] - (thick[1:end-1] + thick[2:end]).*velocity[2:end-1] + thick[2:end].*velocity[3:end])./thick_edgeFace  +
        - 2.0*dt/froude_number_squared*0.5.*(grad_thick[1:end-1] + grad_thick[2:end])
    
    global stress_new[2:end-1] = stress[2:end-1] + dt*froude_upon_mach2*(0.5*thick_edgeFace.*(velocity[3:end]-velocity[1:end-2])/dx  - argand_number*stress[2:end-1]./viscosity[2:end-1]) +
        + 0.5*(dt/dx)^2/mach_number_squared*(stress[1:end-2]./thick[1:end-1] - (1.0./thick[1:end-1] + 1.0./thick[2:end]).*stress[2:end-1] + stress[3:end]./thick[2:end]).*thick_edgeFace +
        - dt^2/mach_number_squared*thick_edgeFace.*(grad_thick[2:end]-grad_thick[1:end-1])/dx +
        - 0.5*dt^2*froude_upon_mach2^2*argand_number*(0.5*thick_edgeFace.*(velocity[3:end]-velocity[1:end-2])/dx  - argand_number*stress[2:end-1]./viscosity[2:end-1])
    
    global velocity=copy(velocity_new)
    global stress= copy(stress_new)
    
    global stress[[1,end]] = [stress[2],2*thick[end]^2-stress[end-1]]
    global velocity[[1,end]] = [(2*U_inflow-velocity[2]), velocity[end-1]]

    #global velocity_new[2:end-1] = velocity[2:end-1]  
    #    + 0.5*(dt/dx)./(thick_edgeFace)./(froude_number_squared).*(stress[3:end]-stress[1:end-2])  
    #    - 0.25*(dt^2/dx)./(thick_edgeFace) * argand_upon_mach2.*(stress[3:end]./viscosity[3:end]-stress[1:end-2]./viscosity[1:end-2])  
    #    + 0.5*(dt/dx)^2/mach_number_squared.*(thick[1:end-1].*velocity[1:end-2] - (thick[1:end-1] + thick[2:end]).*velocity[2:end-1] + thick[2:end].*velocity[3:end])./thick_edgeFace  
    #    - 2*dt/froude_number_squared*0.5.*(grad_thick[1:end-1] + grad_thick[2:end]); 
    

    #global stress_new[2:end-1] = stress[2:end-1] + dt*froude_upon_mach2*(0.5*thick_edgeFace.*(viscosity[2:end-1].*velocity[3:end]-velocity[1:end-2])/dx  - argand_number*stress[2:end-1]) 
    #    + 0.5*(dt/dx)^2/mach_number_squared*(stress[1:end-2]./thick[1:end-1] .- (1 ./thick[1:end-1] .+ 1 ./thick[2:end]).*stress[2:end-1] + stress[3:end]./thick[2:end]).*thick_edgeFace 
    #    - dt^2/mach_number_squared*thick_edgeFace.*(grad_thick[2:end]-grad_thick[1:end-1])/dx 
    #    - 0.5*dt^2*froude_upon_mach2^2*argand_number*(0.5*thick_edgeFace.*viscosity[2:end-1].*(velocity[3:end]-velocity[1:end-2])/dx  - argand_number*stress[2:end-1]);

    #global velocity = copy(velocity_new)
    #global stress = copy(stress_new)
    =#
    global pos = pos + (velocity[2:end-1]-ua)*dt 
    
    
    seismogram[timeCount,2] = pos[end-1];
    seismogram[timeCount,1] = pos[2];

    seismogram[timeCount,4] = velocity[end-1].*x[end]-ua[end];
    seismogram[timeCount,3] = (velocity[2]-ua[1]);
end



if experiment == 3
    data = load("../Data/seimo1.jld")
    tt= (0:size(seismogram)[1]-1)
    p1 = plot(tt*dt*L0/U00/60,seismogram[:,1]*L0*100,legend=false,ylabel="Displacement (cm)",color=:black,linewidth=:2,grid=false,tickfontsize = 10,title="Grounding Line")
   # plot!(data["tt"],data["seismo"][:,1]*L0,color=:dodgerblue,linewidth=:2,linestyle=:dash)
    plot!(xlims=(0,30))
    plot!(ylims=(-1,14))
    annotate!(2.5, 13, "(a)", color=:black)

    p2 = plot(tt*dt*L0/U00/60,seismogram[:,3]*U00/1e3*seconds_in_year,legend=false,ylabel="Velocity (km/a)",xlabel="Time (minutes)",color=:black,linewidth=:2,grid=false,tickfontsize = 10)
    #plot!(data["tt"],data["seismo"][:,3]*data["T"],color=:dodgerblue,linewidth=:2,linestyle=:dash)
    plot!(xlims=(0,30))
    plot!(ylims=(-5,35))
    annotate!(2.5,32, "(c)", color=:black)
    #annotate!(250, 30, text("Ar=0.0005", :dodgerblue,:left,10))
    #annotate!(250, 35, text("Ar=0.00025", :black,:left,10))

    p3 = plot(tt*dt*L0/U00/60,seismogram[:,2]*L0*100,legend=false,color=:black,linewidth=:2,grid=false,tickfontsize = 10,title="Calving Front")
    #plot!(data["tt"],data["seismo"][:,2]*L0,color=:dodgerblue,linewidth=:2,linestyle=:dash)
    plot!(xlims=(0,30))
    plot!(ylims=(-1,14))
    annotate!(2.5,13, "(b)", color=:black)

    p4 = plot(tt*dt*L0/U00/60,seismogram[:,4]*U00/1e3*seconds_in_year,legend=false,xlabel="Time (minutes)",color=:black,linewidth=:2,grid=false,tickfontsize = 10)
    #plot!(data["tt"],data["seismo"][:,4]*data["T"],color=:dodgerblue,linewidth=:2,linestyle=:dash)
    plot!(xlims=(0,30))
    plot!(ylims=(-5,35))
    annotate!(2.5,32, "(d)", color=:black)

    plot(p1,p3,p2,p4,layout=(2,2))
    savefig("../Figures/elastic_test.pdf")
    #JLD2.save("../Data/seimo1.jld","seismo",seismogram,"tt",tt*dt*L0/U00/60,"L0",L0,"T",U00/1e3*seconds_in_year)

end


if experiment<3 && experiment>0
    title = @sprintf("Mach Number=%2.1f",mach_number)
    p1 = plot(x,velocity_array[:,1],legend=false,ylabel=L"\Delta u/u",xticks = [0,0.5,1],yticks = [-0.2,0.2],color=:black,linewidth=:2,grid=false,tickfontsize = 10,xlims=(0,1),ylims=(-0.2,0.2),xlabelfontsize=11,ylabelfontsize=10,framestyle = :box)
    annotate!(1.25,0.3,title,fontsize=:10)
    annotate!(0.05,0.14,text("t=0",:left,10))
    plot!(xticks = ([0,0.5,1], ["0", "0.5", "1"]))
    p2 = plot(x,velocity_array[:,2],legend=false,ylabel=L"\Delta u/u",xticks = [0,0.5,1],yticks = [-0.2,0.2],color=:black,linewidth=:2,grid=false,tickfontsize = 10,xlims=(0,1),ylims=(-0.2,0.2),xlabelfontsize=11,ylabelfontsize=11,framestyle = :box)
    annotate!(0.05,0.14,text("t=0.5",:left,10))
    plot!(xticks = ([0,0.5,1], ["0", "0.5", "1"]))
    p3 = plot(x,velocity_array[:,4],xlabel="Distance",legend=false,ylabel=L"\Delta u/u",xticks = [0,0.5,1],yticks = [-0.2,0.2],color=:black,linewidth=:2,grid=false,tickfontsize = 10,xlims=(0,1),ylims=(-0.2,0.2),xlabelfontsize=11,ylabelfontsize=11,framestyle = :box)
    annotate!(0.05,0.14,text("t=1.5",:left,10))
    plot!(xticks = ([0,0.5,1], ["0", "0.5", "1"]))

    
    p4 = plot(x,stress_array[:,1],legend=false,
        ylabel=L"$\Delta \sigma_{xx}/\sigma_{xx}$",xticks = [0,0.5,1],yticks = [-0.2,0.2],color=:black,linewidth=:2,grid=false,tickfontsize = 10,xlims=(0,1),ylims=(-0.2,0.2),xlabelfontsize=11,ylabelfontsize=11,framestyle = :box)
    annotate!(0.05,0.14,text("t=0",:left,10))
    plot!(xticks = ([0,0.5,1], ["0", "0.5", "1"]))

    p5 = plot(x,stress_array[:,2],
        legend=false,
        ylabel=L"$\Delta \sigma_{xx}/\sigma_{xx}$",
        xticks = [0,0.5,1],yticks = [-0.2,0.2],
        color=:black,linewidth=:2,
        grid=false,tickfontsize = 10,
        xlims=(0,1),ylims=(-0.2,0.2),
        xlabelfontsize=11,
        ylabelfontsize=11,framestyle = :box)
    annotate!(0.05,0.14,text("t=0.5",:left,10))
    plot!(xticks = ([0,0.5,1], ["0", "0.5", "1"]))

    p6 = plot(x,stress_array[:,4],
        legend=false,
        ylabel=L"$\Delta \sigma_{xx}/\sigma_{xx}$",
        xlabel="Distance",
        xticks = [0,0.5,1],
        yticks = [-0.2,0.2],
        color=:black,
        linewidth=:2,
        grid=false,
        tickfontsize = 10,
        xlims=(0,1),
        ylims=(-0.2,0.2),
        xlabelfontsize=11,
        ylabelfontsize=11,framestyle = :box)
    annotate!(0.05,0.14,text("t=1.5",:left,10))
    plot!(xticks = ([0,0.5,1], ["0", "0.5", "1"]))

    plot(p1,p4,p2,p5,p3,p6,layout = (3, 2),titlefont=10,right_margin=5*Plots.mm,top_margin=6*Plots.mm)
    plot_title = @sprintf("../Figures/MachNumber%2.1f.pdf",mach_number)
    savefig(plot_title)
end