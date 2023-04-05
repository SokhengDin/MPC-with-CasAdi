#include "auto_mecanum.hpp"

int main(int argc, char* argv[])
{
    double Lx = 0.165;
    double Ly = 0.18;
    double wheel_radius = 0.05;
    double dt = 0.1;
    double sim_time = 10.0;
    int prediction_horizons = 10;
    int Niter = 0;

    // Boundary
    // Eigen::VectorXd x_min(3);
    // x_min << -3, -3, -1.57;
    // Eigen::VectorXd x_max(3);
    // x_max <<  3,  3,  1.57;
    // Eigen::VectorXd u_min(4);
    // u_min << -10, -10, -10, -10;
    // Eigen::VectorXd u_max(4);
    // u_max << 10, 10, 10, 10;
    std::vector<double> x_min{-3, -3, -1.57};
    std::vector<double> x_max{3, 3, 1.57};
    std::vector<double> u_min{-10, -10, -10, -10};
    std::vector<double> u_max{10, 10, 10, 10};

    // States and Controls
    Eigen::Vector3d current_states(3);
    Eigen::Vector4d current_controls(4);

    Eigen::Vector3d goal_states(3);
    Eigen::Vector4d goal_controls(4);
    goal_states << 2, 2, 0;
    goal_controls << 10, 10, 10, 10;

    // Get solution
    std::vector<double> sol_x;
    std::vector<double> sol_u;

    auto mpc_controller = std::make_shared<AUTO_MECANUM>(AUTO_MECANUM(
        Lx, Ly, wheel_radius,
        dt, sim_time, prediction_horizons
    ));

    // Setup MPC
    mpc_controller->setup_auto_mecanum();

    // Set boundary
    mpc_controller->set_boundary(x_min, x_max, u_min, u_max);

    while (Niter * dt < sim_time)
    {
        mpc_controller->input_trajectory(
            current_states, current_controls,
            goal_states, goal_controls
        );

        sol_x, sol_u = mpc_controller->optimal_solution();

        current_states(sol_x.data());
        current_controls(sol_u.data());

        current_states = current_states + dt * mpc_controller->forward_kinematic(current_states(2),
                                                                                 current_controls(0),
                                                                                 current_controls(1),
                                                                                 current_controls(2),
                                                                                 current_controls(3));

        std::cout << current_states << std::endl;
        Niter = Niter + 1;
    }

}