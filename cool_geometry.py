# Jens Clausen 
# Simulation of water penetration into porous disc

# Imports 
from dolfin import *
from mshr import *
from random import random
import matplotlib.pyplot as plt
import math 

# Define positions 
centre_x = 0.5
centre_y = 0.5
offset_x = 0.0
offset_y = 0.0
outer_radius = 0.5 
inner_radius = 0.2
crit_radius = 0.4
new_radius = 0.05

# Define mesh 
outerCircle = Circle(Point(centre_x,centre_y), outer_radius)
innerCircle = Circle(Point(centre_x - offset_x,centre_y - offset_y), inner_radius)
circle1 = Circle(Point(centre_x - 0.25,centre_y - 0.25), new_radius)

domain = outerCircle - innerCircle - circle1
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
u0 = Constant(0.1)
bc = DirichletBC(V, u0, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("-2", degree=0)      # sink term
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



# Data Analysis 
import numpy
#import vtk

# The source file
#file_name = "/Users/jensbclausen/Desktop/Thesis/1 - Assignments/3 - Sim/poisson000000.vtu"

# Read the source file.
#reader = vtk.vtkXMLUnstructuredGridReader()
#reader.SetFileName(file_name)
#reader.Update()  # Needed because of GetScalarRange
#output = reader.GetOutput()
#potential = output.GetPointData().GetArray("f_22")












# USEFUL PAGES IN TUTORIAL PDF 
# Boundary Conditions - 92 

