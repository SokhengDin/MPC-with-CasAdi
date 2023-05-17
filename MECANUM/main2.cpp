#include "auto_omni.hpp"

int main(int argc, char* argv[])
{
    double r = 0.06;
    double L = 0.22;
    double dt = 0.1;
    double sim_time = 20.0;
    int n_states = 3;
    int n_controls = 4;
    int prediction_horizons = 10;
    int Niter = 0;

    // Boundary
    std::vector<double> x_min{-10, -10, -1.57};
    std::vector<double> x_max{ 10,  10, 1.57};
    std::vector<double> u_min{-10, -10, -10, -10};
    std::vector<double> u_max{10, 10, 10, 10};

    // States and Controls
    Eigen::Vector3d current_states(3);
    Eigen::Vector4d current_controls(4);

    Eigen::Vector3d goal_states(3);
    Eigen::Vector4d goal_controls(4);
    goal_states << 10, 10, 1.57;
    goal_controls << 10, 10, 10, 10;

    auto mpc_controller = std::make_shared<AUTO_OMNI>(AUTO_OMNI(
        r, L, dt, prediction_horizons, n_states, n_controls
    ));

    // Setup MPC
    mpc_controller->setup_mpc();

    // Set boundary
    mpc_controller->set_boundary(x_min, x_max, u_min, u_max);

    // Solution
    std::vector<double> results_all;
    std::vector<double> result_x;
    std::vector<double> result_u;

    while (Niter * dt < sim_time)
    {
        auto start_time = std::chrono::high_resolution_clock::now();
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

        current_states = current_states + dt * mpc_controller->forward_kinematic(
                                                                         result_u[0],
                                                                         result_u[1],
                                                                         result_u[2],
                                                                         result_u[3],
                                                                         result_x[2]);

        // auto x_next = current_states + dt * mpc_controller->forward_kinematic(
        //     result_u[0], result_u[1], result_u[2], result_u[3], result_x[2]
        // );

        // current_states = x_next;

        // std::cout << result_x.size() << std::endl;
        std::cout << "Current states: " << current_states << std::endl;
        // std::cout << results_all << std::endl;
        // std::cout << "Current controls: " << current_controls << std::endl;
        // std::cout << "Optimal solution: " << result_x[0] << " " << result_x[1]  << " " << result_x[2] << std::endl << std::endl;
        std::cout << "Optimal control: " << result_u[0] << " " << result_u[1] << " " << result_u[2] << " " << result_u[3] << std::endl << std::endl;
        // std::cout <<  Niter << std::endl;

        auto stop_time = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time);

        // std::cout << "average calculation time for each iteration [s]: "
        //       << duration.count() / 1e6 << std::endl;

        Niter = Niter + 1;
    }

}
