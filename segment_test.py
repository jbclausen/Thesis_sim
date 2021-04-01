# Jens Clausen 
# Simulation of water penetration into porous disc

# Imports 
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import math 
import numpy as np

# Time-step Constants 
height = 0.00266+0.001 # In time step code, need check that height < 2*radius
dt = 0.01

# Physical Constants 
K = 0.5          # Permeability of the disc 
mu = 0.983E-3    # Viscosity of Fluid
r_p = 1E-5       # Pore Radius 

# Define mesh 
centre_x = 0.0
centre_y = 0.0
radius = 0.0127
source_height = 0.00226
outerCircle = Circle(Point(centre_x,centre_y), radius)
domain = outerCircle
mesh = generate_mesh(domain, 200)

V = FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundary
def boundary(x):
    boundary_pt = False 
    if (x[1] <= (-radius + source_height)):
        boundary_pt = True
    return boundary_pt

def boundary1(x):
    boundary_pt = False 
    #if (x[1] >= (-radius + 0.00366)):
    if (x[1] >= (-radius + height)):
        boundary_pt = True
    return boundary_pt

# Define boundary condition
u0 = Constant(0)
bc = DirichletBC(V, u0, boundary)
u1 = Constant(-4800)
bc1 = DirichletBC(V, u1, boundary1)
bcs = [bc, bc1]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("0", degree=0)
g = Expression("0", degree=0)
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Save solution in ParaView format
file = File("poisson.pvd")
file << u



# Get nodal values 
u_nodal_values = u.vector()
u_coor = V.mesh().coordinates()
u_array = np.array(u_nodal_values)
print("Soln shape", u_array.shape)

# Get coordinates 
coor = mesh.coordinates()
coor = V.mesh().coordinates()

# Computation of grad(u)
grad_u2 = project(grad(u), VectorFunctionSpace(mesh, 'Lagrange', 1))
grad_u = np.array(grad_u2.vector())
print(grad_u)

max_height = -radius
dist_between_nodes = 10*radius / sqrt(u_array.shape[0]) # Estimate distance between nodes 
count = 0

grad_sum = 0
height_sum = 0
grad_sum_1 = 0
grad_sum_2 = 0

max_height = height
# Print values 
if coor.shape[0] == u_array.shape[0]:  # degree 1 element

    #for i in range(len(u_array)):
    #    if (grad_u[2*i] < 1E4 and grad_u[2*i] > -1E4 and grad_u[2*i+1] < 1E4 and grad_u[2*i+1] > -1E4):
    #        grad_sum_1 = grad_sum_1 + grad_u[2*i]
    #        grad_sum_2 = grad_sum_2 + grad_u[2*i+1]

    for i in range(len(u_array)):
        point = (coor[i][0], coor[i][1])

        if (u(point) > -4000 and u(point) < -3500 and coor[i][0] > -0.005 and coor[i][0] < 0.005):
            count = count + 1

            #vel = - K * grad_u[2*i+1] / mu # From Darcy Law 
            vel = - r_p*r_p*grad_u[2*i+1]/(32*(height-source_height))     # From Poisouille FLow 
            print(grad_u[2*i+1])
            #new_h = coor[i][1] + dt*vel 
            new_h = height + dt*vel

            grad_sum = grad_sum + grad_u[2*i+1]
            
            #print('u(%8g,%8g) = %g, grad_u = (%8g, %8g)' % 
                #(coor[i][0], coor[i][1], u(point), grad_u[2*i], grad_u[2*i+1]))
            
            height_sum = height_sum + new_h

            if (new_h > max_height):
                max_height = new_h
                #print('u(%8g,%8g) = %g, grad_u = (%8g, %8g)' % 
                    #(coor[i][0], coor[i][1], u(point), grad_u[2*i], grad_u[2*i+1]))

print("Height", height)
print("Max Height", max_height)
print("Av. Height", height_sum/count)
print("Count", count)
print("Av Grad", grad_sum/count)

print("grad x sum", grad_sum_1)
print("grad y sum", grad_sum_2)


#u_P1 = project(u, V)
#u_nodal_values = u_P1.vector()


# Plot solution
soln = plot(grad(u))
plt.title("Segment Case Simulation")
plt.ylabel("y position (m)")
plt.xlabel("x position (m)")
cbar = plt.colorbar(soln)
cbar.set_label("Pressure Gradient")
plt.show()

# Plot mesh 
#plot(mesh)
#plt.show()

# A sink term of -52500000 gives -p_c at the boundary 

#print('u(%8g,%8g) = %g, grad_u = (%8g, %8g)' % 
        #(coor[i][0], coor[i][1], u_array[i], grad_u[2*i], grad_u[2*i+1]))
