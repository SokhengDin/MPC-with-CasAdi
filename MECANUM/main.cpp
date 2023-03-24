#include "include/mpc_mecanum.hpp"


int main()
{
    double Lx = 0.165;
    double Ly = 0.18;
    double wheel_radius = 0.05;
    double time_step = 0.1;
    double sim_time = 20.0;
    const int prediction_horizon = 20;

    auto mpc_controller = std::make_shared<MPCMECANUM>(MPCMECANUM(
        Lx, Ly, wheel_radius, time_step, sim_time, prediction_horizon
    ));

    mpc_controller->setup_mpc();

    // Set boundary
    Eigen::Vector4d u_min;
    u_min << -10, -10, -10, -10;
    Eigen::Vector4d u_max;
    u_max << 10, 10, 10, 10;
    Eigen::Vector3d x_min;
    x_min << -3, -3, -3.14;
    Eigen::Vector3d x_max;
    x_max << 3, 3, 3.14;

    mpc_controller->set_boundary(u_min, u_max, x_min, x_max);

    // Setup simulation parameters

    int mpciter = 0;

    Eigen::Vector3d current_state(1, 3);
    current_state << 0.0, 0.0, 0.0;
    Eigen::Vector3d target_state(1, 3);
    target_state << 3.0, 3.0, 0.0;
    Eigen::Vector4d current_control(1, 4);
    current_control << 0.0, 0.0, 0.0, 0.0;
    Eigen::Vector4d target_control(1, 4);
    target_control << 10, 10, 10, 10;

    std::vector<double> sol_x;
    std::vector<double> sol_u;

    Eigen::VectorXd forward_kinematic = mpc_controller->kinematic_eigen(current_state, current_control);

    while (mpciter * time_step < sim_time)
    {
        mpc_controller->trajectory_generation(
            current_state, target_state,
            current_control, target_control
        );

        sol_x, sol_u = mpc_controller->get_optimal_solution();

        current_state = sol_x[0];
    }
}