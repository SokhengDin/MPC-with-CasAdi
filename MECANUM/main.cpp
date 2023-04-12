#include "auto_mecanum.hpp"

int main(int argc, char* argv[])
{
    double Lx = 0.165;
    double Ly = 0.18;
    double wheel_radius = 0.05;
    double dt = 0.1;
    double sim_time = 10.0;
    int prediction_horizons = 100;
    int Niter = 0;

    // Boundary
    std::vector<double> x_min{-3, -3, -1.57};
    std::vector<double> x_max{3, 3, 1.57};
    std::vector<double> u_min{-10, -10, -10, -10};
    std::vector<double> u_max{10, 10, 10, 10};

    // States and Controls
    Eigen::Vector3d current_states(3);
    Eigen::Vector3d states(3);
    Eigen::Vector4d current_controls(4);

    Eigen::Vector3d goal_states(3);
    Eigen::Vector4d goal_controls(4);
    goal_states << 2, 2 , 1.57;
    goal_controls << 10, 10, 10, 10;

    auto mpc_controller = std::make_shared<AUTO_MECANUM>(AUTO_MECANUM(
        Lx, Ly, wheel_radius,
        dt, sim_time, prediction_horizons
    ));

    // Setup MPC
    mpc_controller->setup_auto_mecanum();

    // Set boundary
    mpc_controller->set_boundary(x_min, x_max, u_min, u_max);

    // Solution
    std::vector<double> results_all;
    std::vector<double> result_x;
    std::vector<double> result_u;

    while (Niter * dt < sim_time)
    {
        mpc_controller->input_trajectory(
            current_states, current_controls,
            goal_states, goal_controls
        );

        results_all = mpc_controller->optimal_solution();

        result_x.assign(results_all.begin(), results_all.begin()+(prediction_horizons+1)*3);
		result_u.assign(results_all.begin()+(prediction_horizons+1)*3,
                        results_all.begin()+(prediction_horizons+1)*3 + (prediction_horizons)*4);

        current_states << result_x[0], result_x[1], result_x[2];
        current_controls << result_u[0], result_u[1], result_u[2], result_u[3];

        states = current_states + dt * mpc_controller->forward_kinematic(result_x[2],
                                                                         result_u[0],
                                                                         result_u[1],
                                                                         result_u[2],
                                                                         result_u[3]);
        current_states = states;

        // std::cout << result_x.size() << std::endl;
        std::cout << "Current states: " << current_states << std::endl;
        // std::cout << results_all << std::endl;
        // std::cout << "Current controls: " << current_controls << std::endl;
        // std::cout << "Optimal solution: " << result_x[0] << " " << result_x[1]  << " " << result_x[2] << std::endl << std::endl;
        // std::cout << "Optimal control: " << result_u[0] << " " << result_u[1] << " " << result_u[2] << " " << result_u[3] << std::endl << std::endl;
        // std::cout <<  Niter << std::endl;

        Niter = Niter + 1;
    }

}
