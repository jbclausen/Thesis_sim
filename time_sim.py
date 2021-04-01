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

height_changes = []

# Perform Simulation 
time = 0
while time >= 0:

    print("Time step: \t", time, "\tHeight: \t", height)

    # End Simulation once the front is at top of disc 
    if (height > 2*radius - DOLFIN_EPS):
        height = 2*radius
        break
    
    # Define Dirichlet boundary for p = -p_{CAP}
    def boundary(x):
        boundary_pt = False 
        if (x[1] >= (-radius + height)):
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
    f = Expression("x[1] < (- radius + height + DOLFIN_EPS) ? p1 : 0", p1=p1, radius=radius, height=height, degree=1)      # sink term
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

    # Plot and save solution
    soln = plot(u)
    plt.title(("Segment Case Simulation at Time Step Number %d") % time)
    plt.ylabel("y position (m)")
    plt.xlabel("x position (m)")
    cbar = plt.colorbar(soln)
    cbar.set_label("Pressure (Pa)")
    plt.savefig(("./TimeSim/Pics/TimeSim-%d.png" % time)) # save plot as png 
    #plt.show()
    plt.clf() # clear plot 

    count = 0
    height_sum = 0
    if coor.shape[0] == u_array.shape[0]:  # degree 1 element

        for i in range(len(u_array)):
            point = (coor[i][0], coor[i][1])

            if (u(point) > -4000 and u(point) < -3500 and coor[i][0] > -0.005 and coor[i][0] < 0.005):
                count = count + 1

                vel = - r_p*r_p*grad_u[2*i+1]/(32*(height-source_height))     # From Poisouille FLow 

                new_h = height + dt*vel
                height_sum = height_sum + new_h

    if (count == 0):
        temp = height
        height = 2*height - old_height
        old_height = temp
    else: 
        temp = height
        height = height_sum/count 
        old_height = temp 

    height_changes.append(height - old_height)
    print(count, height_changes)

    time = time + 1

print("Simulation Done!")