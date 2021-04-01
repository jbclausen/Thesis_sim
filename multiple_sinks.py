# Jens Clausen 
# Simulation of water penetration into porous disc with multiple sink terms 

# Imports 
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import math 

crit_radius = 0.4

# Define positions 
centre_x = 0.5
centre_y = 0.5
offset_x = 0.0
offset_y = 0.0
outer_radius = 0.5 
inner_radius = 0.2

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
u0 = Constant(0.14)
bc = DirichletBC(V, u0, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

p1 = Expression("-2", degree=0, domain=mesh)
p2 = Expression("-10", degree=0, domain=mesh)

f = Expression("sqrt((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5)) < 0.2 + DOLFIN_EPS ? 0 : (sqrt((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5)) < 0.35 + DOLFIN_EPS ? p2 : 0)",
               p1=p1, p2=p2, degree=1)      # sink term
g = Expression("0", degree=0)       # no crossing Neumann Boundary 
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Save solution in VTK format
file = File("poisson.pvd")
file << u

# Plot solution
soln = plot(u)
plt.colorbar(soln)
plt.show()