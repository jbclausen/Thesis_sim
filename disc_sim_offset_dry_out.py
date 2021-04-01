# Jens Clausen 
# Simulation of water penetration into porous disc

# Imports 
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import math 

# Define positions 
centre_x = 0.0
centre_y = 0.0
offset_x = 0.006
offset_y = 0.006
outer_radius = 0.0127 
inner_radius = 0.003

# Define mesh 
outerCircle = Circle(Point(centre_x,centre_y), outer_radius)
innerCircle = Circle(Point(centre_x - offset_x,centre_y - offset_y), inner_radius)
domain = outerCircle - innerCircle
mesh = generate_mesh(domain, 200)

V = FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    boundary_pt = False 

    # Interior circle   
    if (math.sqrt(math.pow(x[0]-(centre_x-offset_x),2) + math.pow(x[1]-(centre_y-offset_y),2)) <= (inner_radius + DOLFIN_EPS)):
        boundary_pt = True

    return boundary_pt

# Define boundary condition
u0 = Constant(0.37)
bc = DirichletBC(V, u0, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("-0.3526631176470588*63500000", degree=0)      # sink term
g = Expression("0", degree=0)       # no crossing Neumann Boundary 
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
plt.title("Offset Case Simulation")
plt.ylabel("y position (m)")
plt.xlabel("x position (m)")
cbar = plt.colorbar(soln)
cbar.set_label("Pressure (Pa)")
plt.show()

# Plot mesh 
#plot(mesh)
#plt.title("Offset Case Mesh Visualization")
#plt.ylabel("y position (m)")
#plt.xlabel("x position (m)")
#plt.show()

# A sink term of -63500000 gives -p_c at the boundary 




