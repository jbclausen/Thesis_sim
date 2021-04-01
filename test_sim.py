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
offset_x = 0.0
offset_y = 0.0
outer_radius = 0.0127 
inner_radius = 0.003

# Define mesh 
outerCircle = Circle(Point(centre_x,centre_y), outer_radius)
innerCircle = Circle(Point(centre_x - offset_x,centre_y - offset_y), inner_radius)
domain = outerCircle - innerCircle
mesh = generate_mesh(domain, 200)

V = FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundary 
def boundary(x):
    boundary_pt = False 

    # Interior circle   
    if (math.sqrt(math.pow(x[0]-(centre_x-offset_x),2) + math.pow(x[1]-(centre_y-offset_y),2)) <= (inner_radius + DOLFIN_EPS)):
        boundary_pt = True

    # Exterior circle
    #if (math.sqrt(math.pow(x[0]-(centre_x-offset_x),2) + math.pow(x[1]-(centre_y-offset_y),2)) >= (outer_radius - DOLFIN_EPS - 0.001)):    
        #boundary_pt = True
    
    return boundary_pt

def boundary2(x):
    boundary_pt = False 
    
    if (math.sqrt(math.pow(x[0]-(centre_x-offset_x),2) + math.pow(x[1]-(centre_y-offset_y),2)) >= (outer_radius - DOLFIN_EPS - 0.001)):    
        boundary_pt = True
    
    return boundary_pt

# Define boundary condition
u0 = Constant(0)
u1 = Constant(0)
bc = DirichletBC(V, u0, boundary)
bc2 = DirichletBC(V, u1, boundary2)
bcs = [bc, bc2]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("-2.6*3*59952730", degree=0)      # sink term
g = Expression("0", degree=0)          # no crossing Neumann Boundary 
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Save solution in ParaView format
file = File("poisson.pvd")
file << u

# Plot solution
soln = plot(grad(u))
plt.title("1-D Case Simulation")
plt.ylabel("y position (m)")
plt.xlabel("x position (m)")
cbar = plt.colorbar(soln)
cbar.set_label("Pressure (Pa)")
plt.show()

print(soln)

# Plot mesh 
#plot(mesh)
#plt.title("1-D Case Mesh Visualization")
#plt.ylabel("y position (m)")
#plt.xlabel("x position (m)")
#plt.show()

# sink term of -170000000 gives -p_c at the boundary 
# sink term of -59952730 gives -p_c = -4800 at the boundary

# sink of -2.6*3*59952730 for source in interior and exterior gives limit of dry out 

# interesting code: 
# https://fenicsproject.org/qa/13464/verifying-elasticity-benchmark-structure-boundary-conditions/
