#include <dolfin.h>
#include <FenicsSolidMechanics.h>
#include "../forms/p2_forms/etrial.h"
#include <dolfin/la/PETScVector.h>

using namespace dolfin;

// Wave on left end
class Pulse : public Expression
{
public:
  //Constructor
  Pulse(const double& t) : t(t), Expression(2) {}
  // Gaussian Pulse
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 0.0;
    values[1] = -0.01*sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1])*exp(-1000.*(pow(x[0]-0.5,2.)+pow(x[1]-0.5,2.) - t));
  }
  //Current time
  const double& t;
};

// Right Boundary
class RightBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return std::abs(x[0] - 1.0) < DOLFIN_EPS; }
};

// Left Boundary
class LeftBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return x[0] < DOLFIN_EPS; }
};

// Top Boundary
class TopBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return std::abs(x[1] - 1.0) < DOLFIN_EPS; }
};

// Bottom Boundary
class BottomBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return x[1] < DOLFIN_EPS; }
};

int main(int argc, char *argv[])
{
  printf("%s\n%s\n", argv[1], argv[2]);
  Timer timer("Total plasicity solver time");
  set_log_level(ERROR);

  // Create mesh
  UnitSquareMesh mesh(324, 324);

  // Young's modulus and Poisson's ratio
  const double E = 500000.0; //5000000; //150000;  // N/m^2 or Pa
  const double nu = 0.35;
  double cohesion_random, phi_random;
  sscanf(argv[1],"%lf",&cohesion_random);
  sscanf(argv[2],"%lf",&phi_random);
  const double cohesion = cohesion_random;   // 16000 N/m^2 or Pa
  const double dilatancy_angle = phi_random;   // Radians --- this will be 20 degrees
  const double hardening_parameter = 100.0; //0.0001
  const double friction_angle = dilatancy_angle;   // Radians ----- this will be 20 degrees

  // Temporal discretization
  double rho = 1750.0;   //3000   kg/m^3
  double beta = 0.25;
  double dt  = 0.001; // 0.0000001;
  double eta = 0.2;
  Constant _rho(rho);
  Constant _beta(beta);
  Constant _dt(dt);
  Constant _eta(eta);

  // External forces
  Constant zero(0.0, 0.0);

  // Function spaces
  etrial::FunctionSpace V(mesh);

  // Extract elements for stress and tangent
  std::shared_ptr<const FiniteElement> element_t;
  std::shared_ptr<const FiniteElement> element_s;

  etrial::CoefficientSpace_t Vt(mesh);
  element_t = Vt.element();
  etrial::CoefficientSpace_s Vs(mesh);
  element_s = Vs.element();

  // Create boundary condition
  double t = 0.0;
  Pulse U(t);
  LeftBoundary left;      // auto left = std::make_shared<LeftBoundary>();
  TopBoundary top;        // auto top = std::make_shared<TopBoundary>();
  RightBoundary right;    // auto right = std::make_shared<RightBoundary>();
  BottomBoundary bottom;  // auto bottom = std::make_shared<BottomBoundary>();

  // Setting Dirichlet BC
  DirichletBC bcl(V, zero, left);   // Reflective DBC
  DirichletBC bct(V, zero, top);    // Reflective DBC    DirichletBC bct(V3, c, top); 
  DirichletBC bcr(V, zero, right);  // Reflective DBC    DirichletBC bcr(V1, c, right); 
  DirichletBC bcb(V, zero, bottom); // Reflective DBC    DirichletBC bcb(V2, c, bottom);   

// Collect all DBC in one container 'bcs' 
  std::vector<const DirichletBC*> bcs;
  bcs.push_back(&bcl); bcs.push_back(&bcr);   // std::vector<const DirichletBC*> bcs = {{&bcl, &bcr, &bcb, &bct}};
  bcs.push_back(&bcb); bcs.push_back(&bct);

  // Functions
  Function u(V); //auto u = std::make_shared<Function>(V);  
  Function u0(V); // auto u0 = std::make_shared<Function>(V);    
  Function v(V);  // auto v = std::make_shared<Function>(V);   
  Function a(V);
  u0.interpolate(U);
// Set initial conditions and initialise acceleration function
  // u0->vector()->zero();

  // Object of class DP
  const fenicssolid::DruckerPrager DP(E, nu, friction_angle, dilatancy_angle, cohesion, hardening_parameter);

  // Constituive update
  std::shared_ptr<fenicssolid::ConstitutiveUpdate>
    constitutive_update(new fenicssolid::ConstitutiveUpdate(u, *element_s, *Vs.dofmap(), DP));

  // Create forms and attach functions
  fenicssolid::QuadratureFunction tangent(mesh, element_t,
                                          constitutive_update,
                                          constitutive_update->w_tangent());

  etrial::BilinearForm J(V, V);
  J.t = tangent;
  J.dt = _dt; J.beta = _beta; J.rho = _rho; J.eta = _eta;

  etrial::LinearForm L(V);
  L.f = zero; L.g = zero;
  L.dt = _dt; L.beta = _beta; L.rho = _rho; L.u0 = u0; L.u = u; L.eta = _eta; L.v = v; L.a = a;
  fenicssolid::QuadratureFunction stress(mesh, element_s,
                                         constitutive_update->w_stress());
  L.s = stress;

  // Create PlasticityProblem
  fenicssolid::PlasticityProblem nonlinear_problem(J, L, u, tangent,
                                                   stress, bcs, DP);

  // Create nonlinear solver and set parameters
  dolfin::NewtonSolver nonlinear_solver;
  nonlinear_solver.parameters["convergence_criterion"] = "residual";
  nonlinear_solver.parameters["maximum_iterations"]    = 500;
  nonlinear_solver.parameters["relative_tolerance"]    = 1e-6;
  nonlinear_solver.parameters["absolute_tolerance"]    = 1e-15;
  nonlinear_solver.parameters["linear_solver"]         = "gmres";
//  nonlinear_solver.parameters["preconditioner"]        = "hypre_amg";
//  nonlinear_solver.parameters["report"]                = false;

  // File names for output
  char path[200];
  sprintf(path,"/home/crg/Documents/TrialResultsDP_C%lf_PHI%lf/disp.xdmf",cohesion_random, phi_random);
  File file1(path);
  sprintf(path,"/home/crg/Documents/TrialResultsDP_C%lf_PHI%lf/p.xdmf",cohesion_random, phi_random);
  File file2(path);

  // Equivalent plastic strain for visualisation
  CellFunction<double> eps_eq(mesh);

  unsigned int step = 0;
  unsigned int steps = 20;

  while (step < steps)
  {
    // Load parameter
    std::cout << "Time: " << t << std::endl;

    // First step
    if (near(t, 0.0)){
    } else
    
    // Predictions
    *u0.vector() = *u.vector();
    u0.vector()->axpy(dt, *v.vector());
    u0.vector()->axpy((1-2*beta)/2*dt*dt, *a.vector());
    v.vector()->axpy(dt/2, *a.vector());

    // Elasto-plastic calculation
    nonlinear_solver.solve(nonlinear_problem, *u.vector());

    // Corrections
    *a.vector() = *u.vector();
    a.vector()->axpy(-1, *u0.vector());
    *a.vector() *= 1/(beta*dt*dt);
    v.vector()->axpy(dt/2, *a.vector());

    // Update variables
    constitutive_update->update_history();

    // Write output to files
    file1 << u;
    constitutive_update->eps_p_eq().compute_mean(eps_eq);
    file2 << eps_eq;

    t += dt;
    step++;
  }
  cout << "Solution norm: " << u.vector()->norm("l2") << endl;

  timer.stop();
  dolfin::list_timings();
  return 0;
}
