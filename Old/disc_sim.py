# Jens Clausen 
# Simulation of water penetration into porous disc

# Imports 
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import math 


centre_x = 0.5
centre_y = 0.5
offset_x = 0.0
offset_y = 0.0

# Define mesh 
#mesh = UnitSquareMesh(32, 32)
#square = Rectangle(Point(0,0), Point(1.0,1.0))
square = Circle(Point(centre_x,centre_y), 0.5)
#square = UnitSquareMesh(32, 32)
circle = Circle(Point(centre_x - offset_x,centre_y - offset_y), 0.2)
domain = square - circle
mesh = generate_mesh(domain, 200)

V = FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    #return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS
    boundary_pt = False 
    
    # Boundaries 
    #if (x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS): 
        #boundary_pt = True
    
    # Outer circle
    if (math.sqrt(math.pow(x[0]-centre_x,2) + math.pow(x[1]-centre_y,2)) >= (0.5 - DOLFIN_EPS)):
        boundary_pt = True

    # Interior circle   
    if (math.sqrt(math.pow(x[0]-(centre_x-offset_x),2) + math.pow(x[1]-(centre_y-offset_y),2)) <= (0.2 + DOLFIN_EPS)):
        boundary_pt = True

    return boundary_pt

# Define boundary condition
u0 = Constant(0.1)
bc = DirichletBC(V, u0, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
#f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
f = Expression("-2", degree=0)
#g = Expression("sin(5*x[0])", degree=2)
g = Expression("0", degree=0)
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






