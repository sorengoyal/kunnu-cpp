from dolfin import *
import numpy, time
start_time = time.time()
import sys

# Form compiler options
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

def update(u, u0, v0, a0, beta, gamma, dt):
 # Get vectors (references)
 u_vec, u0_vec = u.vector(), u0.vector()
 v0_vec, a0_vec = v0.vector(), a0.vector()
 # Update acceleration and velocity
 a_vec = (1.0/(2.0*beta))*( (u_vec - u0_vec -
 v0_vec*dt)/(0.5*dt*dt) - (1.0-2.0*beta)*a0_vec )
 # v = dt * ((1-gamma)*a0 + gamma*a) + v0
 v_vec = dt*((1.0-gamma)*a0_vec + gamma*a_vec) + v0_vec
 # Update (t(n) <-- t(n+1))
 v0.vector()[:], a0.vector()[:] = v_vec, a_vec
 u0.vector()[:] = u.vector()

# Parameters for the Newmark beta method
beta, gamma = 0.25, 0.5

# Time stepping parameters
dt = 0.0005   #timestep
T = 15.0*dt
t=0.0

# Mesh information
mesh = RectangleMesh(Point(0,0),Point(3,3),300,300)
# Function space definition over the mesh
V = VectorFunctionSpace(mesh, "CG", 2)   #Continuous Galerkin for the displacement field

# Test and Trial function
u1, w = TrialFunction(V), TestFunction(V)

# Material properties
E = float(sys.argv[1])
nu, rho = 0.35, 1800.0
eta = 0.2
lam=E*nu/((1.0+nu)*(1.0-2.0*nu))
mu=E/(2.0*(1.0+nu))

# Initialization of fields (displacement, velocity, acceleration)
u0, v0, a0 = Function(V), Function(V), Function(V)

# Solution at the current time step/test function

def dot_u(u):
  return (gamma/(beta*dt))*(u - u0) - (gamma/beta - 1.0)*v0 - dt*(gamma/(2.0*beta) - 1.0)*a0

def ddot_u(u):
  return (1.0/(beta*dt**2))*(u - u0 - dt*v0) - (1.0/(2.0*beta) - 1.0)*a0

#Initial conditions
f  = Constant((0.0, 0.0))
#zeros = Expression(("0.0","0.0"))
ui  = Expression(("0.0","-40.0*sin(pi*x[0])*sin(pi*x[1])*exp(-1000.*(pow(x[0]-1.5,2)+pow(x[1]-1.5,2) - t))"), t=0.0)
u0 = interpolate(ui, V)
#u0 = interpolate(zeros,V)
#v0 = interpolate(zeros, V)
#a0 = interpolate(zeros, V)

# Stress tensor with damping term
def sigma(u):
    return 2.0*mu*sym(grad(u)) + (lam*tr(grad(u)) + eta*tr(grad(dot_u(u))))*Identity(mesh.geometry().dim())

# Governing equation in the weak form
F = (rho*inner(ddot_u(u1), w) + inner(sigma(u1), sym(grad(w))))*dx - dot(f, w)*dx

# Extract bilinear and linear forms
a, L = lhs(F), rhs(F)

# Applying the Boundary conditions
def boundaries(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, Constant((0.0, 0.0)), boundaries)

# Set up the PDE
u = Function(V)
problem = LinearVariationalProblem(a, L, u, bcs=bc)
solver = LinearVariationalSolver(problem)
solver.parameters["linear_solver"] = "cg"

#Time stepping and solution
file = File("/home/chaitanya/All_Results/RandomRun_Results/Linear_test1/Random_YM" + str(E) + "/RandomYM.pvd")
count=0
while t <= T:
 t += dt
 solver.solve()
 update(u, u0, v0, a0, beta, gamma, dt)
 count+=1
 print count
 if count%3==0 and count>8:
   file << u
