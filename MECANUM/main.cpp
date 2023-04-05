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
		result_u.assign(results_all.begin() + (prediction_horizons+1) * 3, results_all.end());

        current_states(result_x.data());
        current_controls(result_u.data());

        result_x.erase(result_x.begin(), result_x.begin()+3);
		for (int i = 0; i < 3; i++)
		{
			result_x.push_back(*(result_x.end() - 3 + 1));
		}

        result_u.erase(result_u.begin(), result_u.begin()+4);
		for (int i = 0; i < 4; i++)
		{
			result_u.push_back(*(result_u.end() - 4 + 1));
		}

        current_states = current_states + dt * mpc_controller->forward_kinematic(result_x[2],
                                                                                 result_u[0],
                                                                                 result_u[1],
                                                                                 result_u[3],
                                                                                 result_u[4]);
        // std::cout << current_states << std::endl;
        // std::cout << "Optimal Solution: " << result_x[0] << " " << result_x[1]  << " " << result_x[2] << std::endl << std::endl;
        std::cout << "Optimal control: " << result_u[0] << " " << result_u[1] << result_u[2] << " " << result_u[3] << std::endl << std::endl;
        // std::cout <<  Niter << std::endl;
        Niter = Niter + 1;
    }

}