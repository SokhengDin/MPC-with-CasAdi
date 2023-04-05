#include "auto_mecanum.hpp"

AUTO_MECANUM::AUTO_MECANUM(
            double Lx, double Ly, double wheel_radius,
            double dt, double sim_time, int prediction_horizon
        )
{
    Lx_ = Lx;
    Ly_ = Ly;
    wheel_radius_ = wheel_radius;
    dt_ = dt;
    sim_time_ = sim_time;
    prediction_horizon_ = prediction_horizon;
}

AUTO_MECANUM::~AUTO_MECANUM(){}

// Eigen::VectorXd AUTO_MECANUM::forward_kinematic(double theta_, double w1_, double w2_, double w3_, double w4_
Eigen::VectorXd AUTO_MECANUM::forward_kinematic(double theta, double w1, double w2, double w3, double w4)
{
    Eigen::VectorXd inv_vec(4);
    inv_vec << w1, w2, w3, w4;
    Eigen::MatrixXd rot_mat(3,3);
    rot_mat << std::cos(theta), std::sin(theta), 0,
            -std::sin(theta), std::cos(theta), 0,
            0, 0, 1;
    Eigen::MatrixXd J_for(3,4);
    // J_for << 1, -1, -1, 1,
    //          1, 1, -1, -1,
    //          (1/Lx+Ly), (1/Lx+Ly), (1/Lx+Ly), (1/Lx+Ly);
    J_for << 1, 1, 1, 1,
            -1, 1, 1, -1,
            -1/(Lx_+Ly_), 1/(Lx_+Ly_), -1/(Lx_+Ly_), 1/(Lx_+Ly_);

    J_for = (wheel_radius_/4)*J_for;

    Eigen::VectorXd forward_velocity = rot_mat.transpose()*J_for*inv_vec;

    return forward_velocity;
}

// Eigen::VectorXd AUTO_MECANUM::inverse_kinematic(double vx_, double vy_, double vth_)
Eigen::VectorXd AUTO_MECANUM::inverse_kinematic(double theta, double vx, double vy, double vth)
{
    Eigen::VectorXd for_vec(3);
    for_vec << vx, vy, vth;
    Eigen::MatrixXd J_inv(4,3);
    Eigen::MatrixXd rot_mat(3,3);
    rot_mat << std::cos(theta), std::sin(theta), 0,
            -std::sin(theta), std::cos(theta), 0,
            0, 0, 1;
    // J_inv << 1, 1, (Lx+Ly),
    //          -1, 1, (Lx+Ly),
    //          -1, -1, (Lx+Ly),
    //          1, -1, (Lx+Ly);
    J_inv << 1, -1, -(Lx_+Ly_),
            1, 1, (Lx_+Ly_),
            1, 1, -(Lx_+Ly_),
            1, -1, (Lx_+Ly_);

    J_inv = (1/wheel_radius_)*J_inv;

    Eigen::VectorXd inverse_velocity = J_inv*rot_mat*for_vec;

    return inverse_velocity;
}


void AUTO_MECANUM::setup_auto_mecanum()
{
    // States
    MX x = MX::sym("x");
    MX y = MX::sym("y");
    MX theta = MX::sym("theta");
    MX states = MX::vertcat({x, y, theta});
    int n_states = states.size1();
    // Controls
    MX w1 = MX::sym("w1");
    MX w2 = MX::sym("w2");
    MX w3 = MX::sym("w3");
    MX w4 = MX::sym("w4");
    MX controls = MX::vertcat({w1, w2, w3, w4});
    int n_controls = controls.size1();
    // Symbolic Matrix
    MX X = MX::sym("X", n_states, prediction_horizon_+1);
    MX X_ref = MX::sym("X_ref", n_states, prediction_horizon_+1);
    MX U = MX::sym("U", n_controls, prediction_horizon_);
    MX U_ref = MX::sym("U_ref", n_controls, prediction_horizon_);

    // Weight Matrix Diagonalize
    Q_(0,0) = 75;
    Q_(1,1) = 75;
    Q_(2,2) = 90;
    R_(0,0) = 0.01;
    R_(1,1) = 0.01;
    R_(2,2) = 0.01;
    R_(3,3) = 0.01;

    // Right-hand side equation
    /// Rotation Matrix
    rot_mat(0,0) = MX::cos(theta);
    rot_mat(0,1) = MX::sin(theta);
    rot_mat(0,2) = 0;
    rot_mat(1,0) =-MX::sin(theta);
    rot_mat(1,1) = MX::cos(theta);
    rot_mat(1,2) = 0;
    rot_mat(2,0) = 0;
    rot_mat(2,1) = 0;
    rot_mat(2,2) = 1;

    J_for(0,0) = 1;
    J_for(0,1) = 1;
    J_for(0,2) = 1;
    J_for(0,3) = 1;
    J_for(1,0) = -1;
    J_for(1,1) = 1;
    J_for(1,2) = 1;
    J_for(1,3) = -1;
    J_for(2,0) = -1/(Lx_+Ly_);
    J_for(2,1) = 1/(Lx_+Ly_);
    J_for(2,2) = -1/(Lx_+Ly_);
    J_for(2,3) = 1/(Lx_+Ly_);

    RHS.push_back((wheel_radius_/4)*(w1*(MX::sin(theta)+MX::cos(theta))+w2*(-MX::sin(theta)+MX::cos(theta))+w3*(-MX::sin(theta)+MX::cos(theta))+w4*(MX::sin(theta)+MX::cos(theta))));
    RHS.push_back((wheel_radius_/4)*(w1*(MX::sin(theta)-MX::cos(theta))+w2*(MX::sin(theta)+MX::cos(theta))+w3*(MX::sin(theta)+MX::cos(theta))+w4*(MX::sin(theta)-MX::cos(theta))));
    RHS.push_back((wheel_radius_/4)*(-w1+w2+w3-w4)*(1/(Lx_+Ly_)));


    system_kinematic_ = Function("f", {states, controls}, {MX::vertcat(RHS)});

    // std::cout << RHS << std::endl;

    // Create free decision variables
    // std::cout << "Stage 1" << std::endl;
    MX opt_dec = MX::vertcat({MX::reshape(X, -1, 1), MX::reshape(U, -1, 1)});
    MX opt_par = MX::vertcat({MX::reshape(X_ref, -1, 1), MX::reshape(U_ref, -1, 1)});


    // Stage cost
    for (int k = 0; k < prediction_horizon_-1; k++)
    {
        MX st_err = X(all, k) - X_ref(all, k);
        MX con_err = U(all, k) - U_ref(all, k);
        cost_fn = cost_fn + MX::mtimes({st_err.T(), Q_, st_err}) + MX::mtimes({con_err.T(), R_, con_err});
    }

    // std::cout << "Stage 2" << std::endl;

    // Terminal cost
    cost_fn = cost_fn + MX::mtimes({(X(all, prediction_horizon_)-X_ref(all, prediction_horizon_)).T(),
                                    Q_, (X(all, prediction_horizon_)-X_ref(all, prediction_horizon_))});

    // Intialize Constraint
    g = X(all, 0) - X_ref(all, 0);

    // std::cout << "Stage 3" << std::endl;

    // Path Constraint
    for (int k = 0; k < prediction_horizon_-1; k++)
    {
        MX st_next = X(all, k+1);
        // MX st_RK4 = RK4_function((MX)(X(all, k)), (MX)(U(all,k)), system_kinematic_);
        std::vector<MX> input;
        input.push_back(X(all, k));
        input.push_back(U(all, k));
        MX st_euler = system_kinematic_(input).at(0) * dt_ + X(all, k);
        g = MX::vertcat({g, st_next-st_euler});
        // std::cout << st_euler << std::endl;
    }

    // Nonlinear problem
    MXDict nlp_prob = {
        {"f", cost_fn},
        {"x", opt_dec},
        {"p", opt_par},
        {"g", g}
    };

    // Solver
    std::string solver_name = "ipopt";
    Dict nlp_opts;
    nlp_opts["expand"] = true;
    nlp_opts["ipopt.max_iter"] = 1000;
    nlp_opts["ipopt.print_level"] = 0;
    nlp_opts["print_time"] = 0;
    nlp_opts["ipopt.acceptable_tol"] = 1e-8;
    nlp_opts["ipopt.acceptable_obj_change_tol"] = 1e-6;

    solver_ = nlpsol("nlpsol", solver_name, nlp_prob, nlp_opts);

    // std::cout << "Stage 4" << std::endl;
}

// void AUTO_MECANUM::set_boundary(Eigen::Vector3d x_min, Eigen::Vector3d x_max,
//                                 Eigen::Vector4d u_min, Eigen::Vector4d u_max)
// {
//     int num_states = 3;
//     int num_controls = 4;
//     int con_states =  3*(prediction_horizon_+1);
//     int con_controls = con_states + 4*prediction_horizon_;


//     std::cout << "Stage 5" << std::endl;

//     for (int k = 0; k < num_states; k++)
//     {
//         for (int j = 0; j < con_states; j=j+num_states)
//         {
//             lbx_[j] = x_min[k];
//             ubx_[j] = x_max[k];
//         }
//     }
//     for (int k = 0; k < num_controls; k++)
//     {
//         for (int j = con_states; j < con_controls; j=j+num_controls)
//         {
//             lbx_[j] = u_min[k];
//             ubx_[j] = u_max[k];
//         }
//     }

//     std::vector<double> lbx(lbx_.data(), lbx_.data()+lbx_.size());
//     std::vector<double> ubx(ubx_.data(), ubx_.data()+ubx_.size());

//     args_ = {
//         {"lbx", lbx},
//         {"ubx", ubx}
//     };

//     std::cout << "Stage 6" << std::endl;

// }

void AUTO_MECANUM::set_boundary(std::vector<double> x_min, std::vector<double> x_max,
                                std::vector<double> u_min, std::vector<double> u_max)
{
    for (int k = 0; k < prediction_horizon_; k++)
    {
        lbx_.push_back(x_min[0]);
        lbx_.push_back(x_min[1]);
        lbx_.push_back(x_min[2]);

        ubx_.push_back(x_max[0]);
        ubx_.push_back(x_max[1]);
        ubx_.push_back(x_max[2]);
    }

    for (int k = 0; k < prediction_horizon_ -  1; k++)
    {
        lbx_.push_back(u_min[0]);
        lbx_.push_back(u_min[1]);
        lbx_.push_back(u_min[2]);
        lbx_.push_back(u_min[3]);

        ubx_.push_back(u_max[0]);
        ubx_.push_back(u_max[1]);
        ubx_.push_back(u_max[2]);
        ubx_.push_back(u_max[3]);

    }

    // std::cout << "Stage 5" << std::endl;

}

void AUTO_MECANUM::input_trajectory(
        Eigen::Vector3d current_states, Eigen::Vector4d current_controls,
        Eigen::Vector3d goal_states, Eigen::Vector4d goal_controls)
{
    // Params
    Eigen::MatrixXd init_states;
    Eigen::MatrixXd init_controls;
    Eigen::MatrixXd next_trajectories;
    Eigen::MatrixXd next_controls;

    init_states = current_states.replicate(1, prediction_horizon_+1).reshaped();
    init_controls = current_controls.replicate(1, prediction_horizon_).reshaped();
    next_trajectories = goal_states.replicate(1, prediction_horizon_+1).reshaped();
    next_controls = goal_controls.replicate(1, prediction_horizon_).reshaped();

    // std::cout << "Stage 7" << std::endl;

    Eigen::VectorXd variables(init_states.rows()+init_controls.rows());
    Eigen::VectorXd params(next_trajectories.rows()+next_controls.rows());

    variables << init_states,
                init_controls;

    params << next_trajectories,
              next_controls;

    // std::cout << params << std::endl;

    std::vector<double> x0(variables.data(), variables.data()+variables.size());
    std::vector<double> p(params.data(), params.data()+params.size());

    args_ = {
        {"lbg", lbg_},
        {"ubg", ubg_},
        {"x0", x0},
        {"p", p}
    };

    // std::cout << args_["x0"] << std::endl;
}


std::vector<double> AUTO_MECANUM::optimal_solution()
{
    results_ = solver_(args_);

    // std::cout << "Stage 8" << std::endl;

    std::vector<double> results_all(results_.at("x"));

    // std::cout << result_u[0] << std::endl;

    return results_all;
}

