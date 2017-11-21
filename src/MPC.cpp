#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// define desired velocity
#define velocity_des 40

// TODO: Set the timestep length and duration
#define PREDICTION_HORIZON 1 //seconds

size_t N = 10;
double dt = PREDICTION_HORIZON/(double)N; //seconds

// retrieve the number of state variables
#define n_state 6

// define the number of actuators
#define n_actuator 2

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

class FG_eval {
  // the indices of the first vars
  vector<double> start_index;

 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) {
    this->coeffs = coeffs;

    // set the starting indices
    // the order is: x, y, psi, v, cte, e_psi, steering, acceleration
    for (int i = 0; i < n_state; i++) {
      start_index.push_back(i*N);
    }
    for (int i = 0; i < n_actuator; i++) {
      start_index.push_back(n_state*N + i*(N-1));
    }
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    // set the initial state
    for (size_t i = 0; i < n_state; i++) {
      fg[1 + start_index[i]] = vars[start_index[i]];
    }

    // set the constrains
    for (size_t t = 1; t < N; t++) {
      // set the current state
      AD<double> x0 = vars[start_index[0] + t - 1];
      AD<double> y0 = vars[start_index[1] + t - 1];
      AD<double> psi0 = vars[start_index[2] + t - 1];
      AD<double> v0 = vars[start_index[3] + t - 1];
      //AD<double> cte0 = vars[start_index[4] + t - 1];
      AD<double> e_psi0 = vars[start_index[5] + t - 1];
      AD<double> steer0 = vars[start_index[6] + t - 1];
      AD<double> accel0 = vars[start_index[7] + t - 1];

      // set the future state
      AD<double> x1 = vars[start_index[0] + t];
      AD<double> y1 = vars[start_index[1] + t];
      AD<double> psi1 = vars[start_index[2] + t];
      AD<double> v1 = vars[start_index[3] + t];
      AD<double> cte1 = vars[start_index[4] + t];
      AD<double> e_psi1 = vars[start_index[5] + t];

      //calculate the current desired y location - f_x0 and the current desired psi - psi_des
      AD<double> f_x0 = 0;
      AD<double> psi_des = 0;
      // calculate f(x0)
      for (int i = 0; i < coeffs.size(); i++)
        f_x0 += coeffs[i] * CppAD::pow(x0,i);
      // calculate f'(x0)
      for (int i = 1; i < coeffs.size(); i++)
        psi_des += i * coeffs[i] * CppAD::pow(x0,i-1);
      // calculate atan(f'(x0))
      psi_des = CppAD::atan(psi_des);

      // update constraints
      fg[1 + start_index[0] + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + start_index[1] + t] = y1 - (y1 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + start_index[2] + t] = psi1 - (psi0 + (v0/Lf) * steer0 * dt);
      fg[1 + start_index[3] + t] = v1 - (v0 + accel0 * dt);
      fg[1 + start_index[4] + t] = cte1 - ((f_x0-y0) + v0 * CppAD::sin(e_psi0) * dt);
      fg[1 + start_index[5] + t] = e_psi1 - ((e_psi0-psi_des) + (v0/Lf) * steer0 * dt);
    }

    // set the cost function
    fg[0] = 0;
    // punish for position error, angular error, and velocity error
    for (size_t t = 0; t < N; t++) {
      fg[0] += CppAD::pow(vars[start_index[4]+t],2);  //cte^2
      fg[0] += CppAD::pow(vars[start_index[5]+t],2);  //e_psi^2
      fg[0] += CppAD::pow(vars[start_index[3]+t] - velocity_des,2);  //velocity_error^2
    }

    // punish for use of actuators.
    for (size_t t = 0; t < N-1; t++) {
      fg[0] += CppAD::pow(vars[start_index[6] + t], 2); //steering
      fg[0] += CppAD::pow(vars[start_index[7] + t], 2); //acceleration
    }

    // punish for strong change in actuators
    for (size_t t = 0; t < N-2; t++) {
      fg[0] += CppAD::pow(vars[start_index[6]+t+1] - vars[start_index[6]+t],2); //steering diff
      fg[0] += CppAD::pow(vars[start_index[7]+t+1] - vars[start_index[7]+t],2); //acceleration diff
    }

  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {
  // set the starting indices
  // the order is: x, y, psi, v, cte, e_psi, steering, acceleration
  for (int i = 0; i < n_state; i++) {
    start_index.push_back(i*N);
  }
  for (int i = 0; i < n_actuator; i++) {
    start_index.push_back(n_state*N + i*(N-1));
  }
}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = n_state*N + n_actuator*(N-1);
  // TODO: Set the number of constraints
  size_t n_constraints = n_state*N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // set the initial state to vars
  for (i = 0; i < n_state; i++) {
    vars[start_index[i]] = state[i];
  }

  // TODO: Set lower and upper limits for variables.
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  for (i = start_index[n_state]; i < start_index[n_state+1]; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }
  for (i = start_index[n_state+1]; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> actuator_vals;
  for (i = 0; i < n_actuator; i++) {
    int index = start_index[n_state+i];
    actuator_vals.push_back(solution.x[index]);
  }

  // set the predicted trajectory
  ptsx.clear();
  ptsy.clear();
  for (i = 0; i < N; i++) {
    int index_x = start_index[0] + i;
    int index_y = start_index[1] + i;
    ptsx.push_back(solution.x[index_x]);
    ptsy.push_back(solution.x[index_y]);
  }



  return actuator_vals;
}
