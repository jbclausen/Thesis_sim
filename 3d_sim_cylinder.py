from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
from vedo.dolfin import plot as vplot


a = 1.0   # inner radius of iron cylinder
b = 1.2   # outer radius of iron cylinder

# Create mesh and define function space
#mesh = UnitSquareMesh(6, 4)
#mesh = UnitCubeMesh(6, 4, 5)
cylinder = Circle(Point(0, 0), b) - Circle(Point(0, 0), a)

V = FunctionSpace(mesh, 'Lagrange', 1)

# Define boundary conditions
u0 = Expression('0', degree=1)

def u0_boundary(x, on_boundary):
    return x[2] < 0.0 + DOLFIN_EPS

bc = DirichletBC(V, u0, u0_boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6000)
a = inner(nabla_grad(u), nabla_grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Plot solution and mesh
vplot(u)
vplot(mesh)

# Dump solution to file in VTK format
file = File('poisson.pvd')
file << u

# Hold plot
#plt.show()
fig = plt.figure()
ax = plt.axes(projection="3d")

plt.show()