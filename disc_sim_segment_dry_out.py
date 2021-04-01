# Jens Clausen 
# Simulation of water penetration into porous disc

# Imports 
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import math 

# Define mesh 
centre_x = 0.0
centre_y = 0.0
radius = 0.0127
outerCircle = Circle(Point(centre_x,centre_y), radius)
domain = outerCircle
mesh = generate_mesh(domain, 200)

V = FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundary
def boundary(x):

    boundary_pt = False 

    # Interior circle   
    if (x[1] <= (-radius + 0.00266)):
        boundary_pt = True

    return boundary_pt

# Define boundary condition
u0 = Constant(0)
bc = DirichletBC(V, u0, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("-0.3526631176470588*52500000", degree=0)
g = Expression("0", degree=0)
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Save solution in ParaView format
file = File("poisson.pvd")
file << u

# Plot solution
soln = plot(u)
plt.title("Segment Case Simulation")
plt.ylabel("y position (m)")
plt.xlabel("x position (m)")
cbar = plt.colorbar(soln)
cbar.set_label("Pressure (Pa)")
plt.show()

# Plot mesh 
#plot(mesh)
#plt.show()

# A sink term of -52500000 gives -p_c at the boundary 


