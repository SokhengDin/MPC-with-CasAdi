
#include "casadi/casadi.hpp"

using namespace casadi;
using namespace std;

int main(int argc, char *argv[])
{

  // A
  int ncol = 5, nrow = 5;
  vector<casadi_int> colind = {0, 3, 6, 8, 10, 12};
  vector<casadi_int> row = {0, 1, 4, 1, 2, 4, 0, 2, 0, 3, 3, 4};
  vector<double> nz = {19, 12, 12, 21, 12, 12, 21, 16, 21, 5, 21, 18};
  DM A(Sparsity(nrow, ncol, colind, row), nz);

  // Right hand side
  DM b = DM::ones(ncol);

  // Type of linear systems
  enum SymType {UNSYM, SYM, PD};

  // All Linear solvers to be tested
  struct Test {
    string solver;
    SymType type;
  };
  vector<Test> tests;
  tests.push_back({"csparse", UNSYM});
  tests.push_back({"csparsecholesky", PD});
  tests.push_back({"lapacklu", UNSYM});
  tests.push_back({"lapackqr", UNSYM});
  tests.push_back({"ma27", SYM});
  tests.push_back({"mumps", UNSYM});
  tests.push_back({"symbolicqr", UNSYM});
  tests.push_back({"qr", UNSYM});
  tests.push_back({"ldl", SYM});
  tests.push_back({"lsqr", UNSYM});

  // Test all combinations
  for (auto s : {UNSYM, SYM, PD}) {
    DM A_test;
    switch (s) {
    case UNSYM:
      cout << "Unsymmetric linear system" << endl;
      A_test = A;
      break;
    case SYM:
      cout << "Symmetric linear system" << endl;
      A_test = A + A.T();
      break;
    case PD:
      cout << "Positive definite linear system" << endl;
      A_test = mtimes(A.T(), A);
      break;
    }
    for (auto t : tests) {
      if (t.type > s) continue; // Cannot be solved
      if (!Linsol::has_plugin(t.solver)) {
        cout << t.solver << " not available" << endl;
        continue;
      }

      // Solver specific options
      Dict opts;
      if (t.solver == "mumps") {
        opts["symmetric"] = s == SYM || s == PD;
        opts["posdef"] = s == PD;
      }

      // Create a solver instance
      Linsol F("F", t.solver, A_test.sparsity(), opts);

      // Solve
      if (F.sfact(A_test.ptr())) casadi_error("'sfact' failed");
      if (F.nfact(A_test.ptr())) casadi_error("'nfact' failed");
      DM x = densify(b);
      if (F.solve(A_test.ptr(), x.ptr(), x.size2())) casadi_error("'solve' failed");

      // Print the solution
      cout << "solution: " << x << " (" <<  t.solver << ")" << endl;
    }
  }

  return 0;
}