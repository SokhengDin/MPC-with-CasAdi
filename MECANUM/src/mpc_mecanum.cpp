#include "mpc_mecanum.hpp"

MPCMECANUM::MPCMECANUM(
    const double Lx, const double Ly, const double wheel_radius,
    const double time_step, const double sim_time, int prediction_horizon
)
{
    Lx_ = Lx;
    Ly_ = Ly;
    wheel_radius_ = wheel_radius;
    time_step_ = time_step;
    sim_time_ = sim_time;
    prediction_horizon_ = prediction_horizon;
}

MPCMECANUM::~MPCMECANUM(){}

Eigen::VectorXd MPCMECANUM::kinematic_eigen(Eigen::VectorXd states, Eigen::VectorXd controls)
{
    Eigen::MatrixXd rot_mat(3,3);
    rot_mat(0,0) = std::cos(states[2]);
    rot_mat(0,1) = std::sin(states[2]);
    rot_mat(0,2) = 0;
    rot_mat(1,0) = -std::sin(states[2]);
    rot_mat(1,1) = std::cos(states[2]);
    rot_mat(1,2) = 0;
    rot_mat(2,0) = 0;
    rot_mat(2,1) = 0;
    rot_mat(2,2) = 1;

    Eigen::MatrixXd J_for(3,4);
    J_for(0,0) = 1;
    J_for(0,1) = -1;
    J_for(0,2) = -1;
    J_for(0,3) = 1;
    J_for(1,0) = 1;
    J_for(1,1) = 1;
    J_for(1,2) = -1;
    J_for(1,3) = -1;
    J_for(2,0) = 1/(Lx_+Ly_);
    J_for(2,1) = 1/(Lx_+Ly_);
    J_for(2,2) = 1/(Lx_+Ly_);
    J_for(2,3) = 1/(Lx_+Ly_);
    J_for = (wheel_radius_/4)*J_for;

    Eigen::MatrixXd rhs;

    rhs = rot_mat*J_for*controls;

    return rhs;
}


void MPCMECANUM::setup_mpc()
{
    // State
    ca::MX x = ca::MX::sym("x");
    ca::MX y = ca::MX::sym("y");
    ca::MX theta = ca::MX::sym("theta");
    ca::MX states = ca::MX::vertcat({x, y, theta});

    // Control
    ca::MX w1 = ca::MX::sym("w1");
    ca::MX w2 = ca::MX::sym("w2");
    ca::MX w3 = ca::MX::sym("w3");
    ca::MX w4 = ca::MX::sym("w4");
    ca::MX controls = ca::MX::vertcat({w1, w2, w3, w4});

    // Symbolic Matrix
    ca::MX X = ca::MX::sym("X", num_states_, prediction_horizon_+1);
    ca::MX X_ref = ca::MX::sym("X_ref", num_states_, prediction_horizon_+1);
    ca::MX U = ca::MX::sym("U", num_controls_, prediction_horizon_);
    ca::MX U_ref = ca::MX::sym("U_ref", num_controls_, prediction_horizon_);

    // Cost function and constraint
    ca::MX opt_dec = ca::MX::vertcat({ca::MX::reshape(X, -1, 1), ca::MX::reshape(X_ref, -1, 1)});
    ca::MX opt_par = ca::MX::vertcat({ca::MX::reshape(U, -1, 1), ca::MX::reshape(U_ref, -1, 1)});

    g_.push_back(X(all, 0)-X_ref(all, 0));

    // Weighted Matrix
    ca::DM Q = ca::DM::zeros(3,3);
    Q(0,0) = 75;
    Q(1,1) = 75;
    Q(2,2) = 90;
    ca::DM R = ca::DM::zeros(4,4);
    R(0,0) = 0.01;
    R(1,1) = 0.01;
    R(2,2) = 0.01;
    R(3,3) = 0.01;

    // Stage cost
    for (int k = 0; k < prediction_horizon_; k++)
    {
        ca::MX st_err = X(all, k) - X_ref(all, k);
        ca::MX con_err = U(all, k) - U_ref(all, k);
        cost_fn_ = cost_fn_ + ca::MX::mtimes({st_err.T(), Q, st_err}) + ca::MX::mtimes({con_err.T(), R, con_err});
    }

    // Terminal cost
    cost_fn_ = cost_fn_ + ca::MX::mtimes({(X(all, prediction_horizon_+1)-X_ref(all, prediction_horizon_+1)).T(), Q, (X(all, prediction_horizon_+1)-X_ref(all, prediction_horizon_+1))});


    // Euler discretization
    for (int k = 0; k < prediction_horizon_; k++)
    {
        ca::MX st_next = X(all, k+1);
        ca::MX st_next_euler = X(all, k) + time_step_ * get_kinematic_model(ca::MX(X(all, k)), ca::MX(U(all, k)));
        g_.push_back(st_next-st_next_euler);
    }

    ca::MXDict nlp_prob = {
        {"f", cost_fn_},
        {"x", opt_dec},
        {"p", opt_par},
        {"g", ca::MX::vertcat(g_)}
    };

    std::string solver_name = "ipopt";
    ca::Dict nlp_opts;
    nlp_opts["expand"] = true;
    nlp_opts["ipopt.max_iter"] = 1000;
    nlp_opts["ipopt.print_level"] = 0;
    nlp_opts["print_time"] = 0;
    nlp_opts["ipopt.acceptable_tol"] = 1e-8;
    nlp_opts["ipopt.acceptable_obj_change_tol"] = 1e-6;

    solver_ = ca::nlpsol("nlpsol", solver_name, nlp_prob, nlp_opts);
}

void MPCMECANUM::set_boundary(Eigen::VectorXd x_min, Eigen::VectorXd x_max,
                              Eigen::VectorXd u_min, Eigen::VectorXd u_max)
{

    // ---- Setup Constraint for States ---- //
    int size_states = num_states_*(prediction_horizon_+1);
    int size_controls = num_states_*(prediction_horizon_+1)+num_controls_*prediction_horizon_;

    std::vector<double> lbx;
    std::vector<double> ubx;

    for (int i = 0; i < num_states_; i++)
    {
        for (int j = 0; j < size_states; j=j+num_states_)
        {
            lbx[j] = x_min[i];
            ubx[j] = x_max[i];
        }
    }

    for (int i = 0; i < num_controls_; i++)
    {
        for (int j = 0; j < size_controls; j=j+num_controls_)
        {
            lbx[j] = u_min[i];
            ubx[j] = u_max[i];
        }
    }

    ca::DM lbg = ca::DM::zeros(num_states_*(prediction_horizon_+1), 1);
    ca::DM ubg = ca::DM::zeros(num_states_*(prediction_horizon_+1), 1);

    args_ = {
        {"lbx", lbx},
        {"ubx", ubx},
        {"lbg", lbg},
        {"ubg", ubg}
    };
}

void MPCMECANUM::trajectory_generation(Eigen::VectorXd init_states, Eigen::VectorXd target_states,
                                       Eigen::VectorXd init_controls, Eigen::VectorXd target_controls)
{
    // Eigen::MatrixXd next_trajectories(init_states.rows(), init_states.cols()+target_states.cols());
    // Eigen::MatrixXd next_controls(init_controls.rows(), init_controls.cols()+target_controls.cols());

    // next_trajectories = next_trajectories.reshaped(-1, 1);
    // next_controls = next_controls.reshaped(-1, 1);

    std::vector<Eigen::VectorXd> next_trajectories;
    std::vector<Eigen::VectorXd> next_controls;

    next_trajectories[0] = init_states;
    next_controls[0] = init_controls;

    for (int k = 1; k < prediction_horizon_; k++)
    {
        next_trajectories.push_back(target_states);
        next_controls.push_back(target_controls);
    }

    args_["x0"] = next_trajectories;
    args_["p"] = next_controls;

}

std::vector<double> MPCMECANUM::get_optimal_solution()
{
    opt_result_ = solver_(args_);

    std::vector<double> sol_all(opt_result_.at("x"));
    std::vector<double> sol_x;
    std::vector<double> sol_u;
    Eigen::MatrixXd sol_x_matrix;
    Eigen::MatrixXd sol_u_matrix;

    sol_x.assign(sol_all.begin(), sol_all.begin() + prediction_horizon_ * num_states_);
    sol_x_matrix = Eigen::MatrixXd::Map(sol_x.data(), num_states_, prediction_horizon_);
    sol_u.assign(sol_all.begin() + prediction_horizon_ * num_states_,
                 sol_all.begin() + prediction_horizon_ * num_states_ + (prediction_horizon_ - 1)*num_controls_);
    sol_u_matrix = Eigen::MatrixXd::Map(sol_u.data(), num_controls_, prediction_horizon_ - 1);

    return sol_all;
}