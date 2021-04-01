# Jens Clausen 
# Simulation of water penetration into porous disc
# R/R_0 = 2 in order to try to reproduce results from paper

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
inner_radius = outer_radius / 3.5

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

    return boundary_pt

# Define boundary condition
u0 = Constant(0)
bc = DirichletBC(V, u0, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("1.25*-170000000*0.3526631176470588", degree=0)      # sink term
g = Expression("0", degree=0)          # no crossing Neumann Boundary 
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds
0
# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Save solution in ParaView format
file = File("poisson.pvd")
file << u

# Plot solution
soln = plot(u)
plt.title("1-D Case Simulation")
plt.ylabel("y position (m)")
plt.xlabel("x position (m)")
cbar = plt.colorbar(soln)
cbar.set_label("Pressure (Pa)")
plt.show()

# Plot mesh 
#plot(mesh)
#plt.title("1-D Case Mesh Visualization")
#plt.ylabel("y position (m)")
#plt.xlabel("x position (m)")
#plt.show()



# sink term of -59952730 gives -p_c = -4800 at the boundary
# sink term of -3.5*-170000000*0.3526631176470588 gives -p_c = -4800 at the boundary R/R_0 = 1.9 
# sink term of -5.3*-170000000*0.3526631176470588 gives -p_c = -4800 at the boundary R/R_0 = 1.65
# sink term of -6.8*-170000000*0.3526631176470588 gives -p_c = -4800 at the boundary R/R_0 = 1.54
# sink term of -59952730 gives -p_c = -4800 at the boundary






# sink term of -170000000 gives -p_c at the boundary for R/R_0 = 4.23 (me/mec = 0.0287)
# sink term of 3.5*-170000000 gives -p_c at the boundary for R/R_0 = 1.9 (me/mec = 0.5)
# sink term of 5.5*-170000000 gives -p_c at the boundary for R/R_0 = 1.65 (me/mec = 1)
# sink term of 7*-170000000 gives -p_c at the boundary for R/R_0 = 1.54 (me/mec = 1.5)
# sink term of 11.5*-170000000 gives -p_c at the boundary for R/R_0 = 1.38 (me/mec = 3)


