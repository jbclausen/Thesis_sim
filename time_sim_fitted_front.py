# time_sim.py 
# Jens Clausen 
# Time Dependent Simulation of water penetration into porous disc

# Imports 
from dolfin import *
from mshr import *
import scipy.io
import matplotlib.pyplot as plt
import math 
import numpy as np
from scipy.optimize import curve_fit

# Curve to fit to the front 
def f_curve(x,a,b):
    return a + b*pow(x,2) #+ c*pow(x,4)

# Time step constant
dt = 1

# Physical Constants 
height = 0.00416 # In time step code, need check that height < 2*radius
source_height = 0.00266
mu = 0.983E-3    # Viscosity of Fluid
r_p = 1E-5       # Pore Radius 
old_height = height

# Define positions 
centre_x = 0.0
centre_y = 0.0
radius = 0.0127 

# Define mesh 
outerCircle = Circle(Point(centre_x,centre_y), radius)
domain = outerCircle 
mesh = generate_mesh(domain, 200)
V = FunctionSpace(mesh, "Lagrange", 1)

# Perform Simulation 
heights = []
time = 0
c_1 = -radius + height 
c_2 = 0
c_3 = 0
sim_Paused = False # Change to True to pick up from previous sim
while time >= 0:

    if (sim_Paused):
        c_1 = 0.0009070081616048515
        c_2 = -20.080441
        time = 91
        sim_Paused = False

    print("Time step: \t", time, "\tConstants:\t", c_1, c_2)
    height = c_1 
    heights.append(height)

    # End Simulation once the front is at top of disc 
    if (c_1 > radius - DOLFIN_EPS):
        break
    
    # Define Dirichlet boundary for p = -p_{CAP}
    def boundary(x):
        boundary_pt = False 
        
        x_pos = x[0]
        y_pos = f_curve(x_pos,c_1,c_2)
        
        #if (x[1] >= (-radius + height) or x[1] >= y_pos):
        if (x[1] >= y_pos):
            boundary_pt = True
        return boundary_pt

    # Define Dirichlet boundary for p = 0
    def boundary2(x):
        boundary_pt = False
        if (x[1] <= (-radius + source_height)):    
            boundary_pt = True
        return boundary_pt

    # Define boundary condition
    u0 = Constant(-4800)
    u1 = Constant(0)
    bc = DirichletBC(V, u0, boundary)
    bc2 = DirichletBC(V, u1, boundary2)
    bcs = [bc, bc2]
    
    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    p1 = Expression("-18514813", degree=0, domain=mesh)
    #f = Expression("x[1] < (- radius + height + DOLFIN_EPS) ? p1 : 0", p1=p1, radius=radius, height=height, degree=1)      # sink term
    f = Expression("0", degree=0)      # sink term
    g = Expression("0", degree=0)      # no crossing Neumann Boundary 
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx + g*v*ds
    
    # Compute solution
    u = Function(V)
    solve(a == L, u, bcs)
    
    # Get nodal values 
    u_nodal_values = u.vector()
    u_coor = V.mesh().coordinates()
    u_array = np.array(u_nodal_values)
    
    # Get coordinates 
    coor = mesh.coordinates()
    coor = V.mesh().coordinates()

    # Compute grad(u)
    grad_u2 = project(grad(u), VectorFunctionSpace(mesh, 'Lagrange', 1))
    grad_u = np.array(grad_u2.vector())


    count = 0
    dy_sum = 0 
    points = []
    if coor.shape[0] == u_array.shape[0]:  # degree 1 element

        for i in range(len(u_array)):
            point = (coor[i][0], coor[i][1])
            y_pred = f_curve(coor[i][0],c_1,c_2)
            #if (u(point) > -4000 and u(point) < -3500 and coor[i][0] > -0.005 and coor[i][0] < 0.005 and coor[i][1] > -radius and coor[i][1] < radius):
            #if (u(point) > -4000 and u(point) < -3500 and coor[i][0] > -3/4*radius and coor[i][0] < 3/4*radius and coor[i][1] > -radius and coor[i][1] < radius):            
            if (u(point) > -4000 and u(point) < -3500 and coor[i][0] > -radius and coor[i][0] < radius and coor[i][1] > -radius and coor[i][1] < radius):
                count = count + 1
                dp = grad_u[2*i+1]
                if (dp > 0):
                    dp = -dp
                vel = - r_p*r_p*dp/(32*((y_pred+radius)-source_height))     # From Poisouille FLow 
                dy = dt*vel 
                dy_sum = dy_sum + dy

                new_point = (coor[i][0], coor[i][1] + dy)
                points.append(new_point)
 
    print(count)

    points = np.array(points)
    
    # Fit quadratic function to the front 
    popt, pcov = curve_fit(f_curve, points[:,0], points[:,1])
    c_1=popt[0].round(6) 
    c_2=popt[1].round(6)
    # Plot Front
    circle1 = plt.Circle((0, 0), radius, color='b', fill=False)
    n = points.shape[0]
    x_plot=np.linspace(-radius,radius, num=n, dtype=float)
    eqn = f'y={c_1}+{c_2}x^2'

    # Plot solution 
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,8), gridspec_kw={'width_ratios': [1, 1.255]})
    soln = plot(u)
    ax2.set_title(("Segment Case Simulation at Time Step Number %d") % time)
    ax2.set_ylabel("y position (m)")
    ax2.set_xlabel("x position (m)")
    cbar = plt.colorbar(soln)
    cbar.set_label("Pressure (Pa)")
    ax2.set_aspect('equal', adjustable='box') # equal axes 

    ax1.set_title("Position of Liquid Front")
    ax1.set_ylabel("y position (m)")
    ax1.set_xlabel("x position (m)")
    ax1.set_xlim(-1.05*radius, 1.05*radius)
    ax1.set_ylim(-1.05*radius, 1.05*radius)
    ax1.plot(x_plot, f_curve(x_plot, *popt), 'r-')
    ax1.add_patch(circle1)
    ax1.set_aspect('equal', adjustable='box') # equal axes 

    plt.savefig(("./TimeSim/Pics/TimeSim-%d.png" % time)) # save plot as png 
    #plt.show()
    plt.clf() # clear plot 
    plt.close() 

    c_1 = height + (dy_sum/count)

    time = time + 1

print(heights)
print("Simulation Done!")



#    # Plot and save solution
#    soln = plot(u)
#    plt.title(("Segment Case Simulation at Time Step Number %d") % time)
#    plt.ylabel("y position (m)")
#    plt.xlabel("x position (m)")
#    cbar = plt.colorbar(soln)
#    cbar.set_label("Pressure (Pa)")
#    plt.gca().set_aspect('equal', adjustable='box') # equal axes 
#    plt.savefig(("./TimeSim/Pics/TimeSim-%d.png" % time)) # save plot as png 
#    #plt.show()
#    plt.clf() # clear plot 
#
#    plt.scatter(points[:,0], points[:,1])
#    plt.xlim(-1.05*radius, 1.05*radius)
#    plt.ylim(-1.05*radius, 1.05*radius)
#    plt.plot(x_plot, f_curve(x_plot, *popt), 'r-')
#    #plt.gca().add_patch(circle1)
#    plt.gca().set_aspect('equal', adjustable='box')
#    plt.title(("Points Assessed in Time Sim at Time Step Number %d") % time)
#    plt.ylabel("y position (m)")
#    plt.xlabel("x position (m)")
#    plt.savefig(("./TimeSim/PointsPics/TimeSim-%d.png" % time))
#    plt.clf()