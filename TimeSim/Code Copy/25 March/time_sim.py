# Jens Clausen 
# Simulation of water penetration into porous disc

# Imports 
from dolfin import *
from mshr import *
import scipy.io
import matplotlib.pyplot as plt
import math 
import numpy as np

# Time step constants 
height = 0.00266 # In time step code, need check that height < 2*radius
source_height = 0.001
dt = 0.001
K = 0.5
mu = 0.983E-3
last_step = False
total_time = 200*dt # number of seconds we want to simulate 

# Define positions 
centre_x = 0.0
centre_y = 0.0
radius = 0.0127 

# Define mesh 
outerCircle = Circle(Point(centre_x,centre_y), radius)
domain = outerCircle 
mesh = generate_mesh(domain, 200)

V = FunctionSpace(mesh, "Lagrange", 1)

for i in range(int(total_time/dt)): 
    
    print("Time step: \t", i, "\t Height: \t", height)

    # Perform Simulation 
    if (last_step):
        break 

    if (height >= 2*radius):
        height = 2*radius
        last_step = True
    
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

    # Computation of grad(u)
    grad_u2 = project(grad(u), VectorFunctionSpace(mesh, 'Lagrange', 1))
    grad_u = np.array(grad_u2.vector())

    max_height = -radius
    count = 0

    grad_sum = 0
    height_sum = 0

    max_height = height

    # Plot solution
    soln = plot(u)
    plt.title(("Segment Case Simulation at Time Step Number %d") % i)
    plt.ylabel("y position (m)")
    plt.xlabel("x position (m)")
    cbar = plt.colorbar(soln)
    cbar.set_label("Pressure (Pa)")
    plt.savefig(("./TimeSim/Pics/TimeSim-%d.png" % i))
    #plt.show()
    plt.clf()

    if coor.shape[0] == u_array.shape[0]:  # degree 1 element
        print(len(u_array))
        for i in range(len(u_array)):
            point = (coor[i][0], coor[i][1])

            if (u(point) > -3000 and u(point) < -2500 and grad_u[2*i + 1] > 0 and grad_u[2*i + 1] < 1E-1):
                
                count = count + 1

                vel = K * grad_u[2*i+1] / mu # From Darcy Law

                new_h = height + dt*vel

                grad_sum = grad_sum + grad_u[2*i+1]

                height_sum = height_sum + new_h

                if (new_h > max_height):
                    max_height = new_h

    av_height = height_sum/count 
    height = av_height


print("Done")



# https://fenicsproject.org/qa/1460/numpy-arrays-from-fenics-data/


