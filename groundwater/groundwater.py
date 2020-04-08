from firedrake import *

m  = 20
Ly = 0.85
dy = Ly/m
mesh = IntervalMesh(m, 0 ,Ly)

# Time definitions
t   = 0.0
end = 150.0
Ntm = 75
dtmeas = end/Ntm
tmeas = dtmeas

# Define Function space on our mesh.
# Initially we will use a continuous linear Lagrange basis
# Try other order, 1, 2, 3
V = FunctionSpace(mesh, "CG", 1)

# Define timestep value
# Dt = 0.00001
# dt = Constant(Dt)

# Define Crank Nicholson parameter theta not zero does not work yet: need other solver?
theta = 0.0

# Define groundwater constants
mpor  = 0.3
sigma = 0.8
Lc    = 0.05
kperm = 1.0e-8
w     = 0.1
R     = 0.000125
nu    = 1.0e-6
g     = 9.81
alpha = kperm/( nu*mpor*sigma )
gam   = Lc/( mpor*sigma )
fac2  = sqrt(g)/( mpor*sigma )
# 
# ncase = 0 Dirichlet bc, ncase = 1 overflow groundwater into canal section with weir equation:
nncase = 1

# Initial condition
h_prev = Function(V)
h_prev.interpolate(Expression("0.00"))

# Create storage for paraview
outfile = File("./Results/groundwater.pvd")

# Write IC to file for paraview
outfile.write(h_prev , t = t )

# Define trial and test functions on this function space
# h will be the equivalent to h^n+1 in our timestepping scheme

h = Function(V)
phi = TestFunction(V)

# Provide intial guess to non linear solve
h.assign(h_prev)

CFL = 0.3
dt = CFL*0.5*dy*dy  # Based on FD estimate; note that dt must be defined before flux, etc 
# Don't understand why. Does this mean that for variable time step one need to redefine F again and again?

def flux ( h , phi , R ):
  return ( alpha * g * h * dot ( grad (h) , grad (phi) ) - (R * phi )/ ( mpor * sigma ) )

if nncase == 0:
  F = ( (h-h_prev)*phi/dt  + theta * flux ( h , phi , R ) + (1-theta)* flux ( h_prev, phi, R) ) *dx
  # Boundary conditions: Condition at Ly satisfied weakly
  bc1 = DirichletBC(V, 0.07, 1)
  h_problem = NonlinearVariationalProblem( F , h , bcs = bc1)
elif nncase == 1:
  if theta == -1.0:
    aa = (h*phi/dt)*dx+(gam*phi*h/dt)*ds(1)
    L2 = ( h_prev*phi/dt - flux ( h_prev, phi, R) ) *dx 
    L = L2+( gam*phi*h_prev/dt-phi*fac2*max_value(2.0*h_prev/3.0,0.0)*sqrt(max_value(2.0*h_prev/3.0,0.0)) )*ds(1)
  elif theta >= 0.0:
    F = ( (h-h_prev)*phi/dt  + theta * flux ( h , phi , R ) + (1-theta)* flux ( h_prev, phi, R) ) *dx
    # Add boundary contributions at y = 0: HERE HERE; does ds1 set it to the contribution at y=0????
    # A  F2 = ( gam*phi*(h-h_prev)/dt +theta*phi*fac2*max_value(2.0*h/3.0,0.0) )*ds(1)
    F2 = ( gam*phi*(h-h_prev)/dt+theta*phi*fac2*max_value(2.0*h/3.0,0.0)*sqrt(max_value(2.0*h/3.0,0.0))+(1-theta)*phi*fac2*max_value(2.0*h_prev/3.0,0.0)*sqrt(max_value(2.0*h_prev/3.0,0.0)) )*ds(1)
    h_problem = NonlinearVariationalProblem( F+F2 , h )
    h_solver = NonlinearVariationalSolver(h_problem)
                
while (t < end):
  # First we increase time
  t += dt
        
  # Print to console the current time
        
  # Use the solver and then update values for next timestep
  if theta == -1.0:
    solve(aa == L, h)
  elif theta >= 0.0:
    h_solver.solve()
                
  h_prev.assign(h)
                
  # Write output to file for paraview visualisation
  if t>tmeas:
    print('Time is %g',t)
    tmeas = tmeas+dtmeas
    outfile.write(h_prev , t = t )
# End while loop










































































































