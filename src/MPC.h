#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {

  // the indices of the first vars
  vector<double> start_index;

 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

  // store the last predicted trajectory
  vector<double> ptsx;
  vector<double> ptsy;
};

#endif /* MPC_H */
