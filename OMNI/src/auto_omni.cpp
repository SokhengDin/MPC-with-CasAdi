#include "auto_omni.hpp"

AUTO_OMNI::AUTO_OMNI(
    double L, double r, double dt,
    int prediction_horizons, int n_states, int n_controls
)
{
    L_ = L;
    r_ = r;
    dt_ = dt;
    prediction_horizons_ = prediction_horizons;
    n_states_ = n_states;
    n_controls_ = n_controls;
    slice_states = Slice(0, n_states_);
    slice_controls = Slice(0, n_controls_);
    lbg_ = DM::zeros((n_states_*(prediction_horizons_+1)), 1);
    ubg_ = DM::zeros((n_states_*(prediction_horizons_+1)), 1);
    a1_ = M_PI/4;
    a2_ = 3*M_PI/4;
    a3_ = 5*M_PI/4;
    a4_ = 7*M_PI/4;

    Q_(0, 0) = 75;
    Q_(1, 1) = 75;
    Q_(2, 2) = 90;
    R_(0, 0) = 1;
    R_(1, 1) = 1;
    R_(2, 2) = 1;
    R_(3, 3) = 1;
}

AUTO_OMNI::~AUTO_OMNI(){}

void AUTO_OMNI::setup_mpc(){
    // States
    MX x = MX::sym("x");
    MX y = MX::sym("y");
    MX yaw = MX::sym("yaw");
    MX states = MX::vertcat({x, y, yaw});
    // Controls
    MX u1 = MX::sym("u1");
    MX u2 = MX::sym("u2");
    MX u3 = MX::sym("u3");
    MX u4 = MX::sym("u4");
    MX controls = MX::vertcat({u1, u2, u3, u4});
    // Symbolics
    MX X = MX::sym("X", n_states_, prediction_horizons_+1);
    MX X_ref = MX::sym("X_ref", n_states_, prediction_horizons_+1);
    MX U = MX::sym("U", n_controls_, prediction_horizons_);
    MX U_ref = MX::sym("U_ref", n_controls_, prediction_horizons_);

    // Forward matrix
    J_for(0, 0) =  MX::sin(yaw+a1_);
    J_for(0, 1) = -MX::sin(yaw+a2_);
    J_for(0, 2) =  MX::sin(yaw+a3_);
    J_for(0, 3) = -MX::sin(yaw+a4_);
    J_for(1, 0) =  MX::cos(yaw+a1_);
    J_for(1, 1) = -MX::cos(yaw+a2_);
    J_for(1, 2) =  MX::cos(yaw+a3_);
    J_for(1, 3) = -MX::cos(yaw+a4_);
    J_for(2, 0) = 1/(2*L_); 
    J_for(2, 1) = 1/(2*L_); 
    J_for(2, 2) = 1/(2*L_); 
    J_for(2, 3) = 1/(2*L_);


    J_for = (r_/2)*J_for;

    rhs_ = MX::mtimes({J_for, controls});

    system_kinematic_ = Function("f", {states, controls}, {rhs_});

    MX opt_var = MX::vertcat({
        MX::reshape(X, -1, 1),
        MX::reshape(U, -1, 1)
    });

    MX opt_dec = MX::vertcat({
        MX::reshape(X_ref, -1, 1),
        MX::reshape(U_ref, -1, 1)
    });


    g_ = X(slice_states, 0) - X_ref(slice_states, 0);
    g_ = MX::reshape(g_, -1, 1);

    for (int k = 0 ; k < prediction_horizons_; k++)
    {
        MX st_err = X(slice_states, k) - X_ref(slice_states, k);
        MX con_err = U(slice_controls, k) - U_ref(slice_controls, k);
        cost_fn_ = cost_fn_ + MX::mtimes({st_err.T(), Q_, st_err}) + MX::mtimes({con_err.T(), R_, con_err});
    }



    cost_fn_ = cost_fn_ + MX::mtimes({(X(slice_states, prediction_horizons_)-X_ref(slice_states, prediction_horizons_)).T(),
                                    Q_,
                                    (X(slice_states, prediction_horizons_)-X_ref(slice_states, prediction_horizons_))});
    
    for (int k = 0; k < prediction_horizons_; k++)
    {
        MX st_next = X(slice_states, k+1);
        std::vector<MX> input;
        input.push_back(X(slice_states, k));
        input.push_back(U(slice_controls, k));
        // MX st_RK4 = RK4_function((MX)(X(slice_states, k)), (MX)(U(slice_controls,k)), system_kinematic_);
        MX st_euler = X(slice_states, k) + system_kinematic_(input).at(0) * dt_;
        g_ = MX::vertcat({g_, st_next-st_euler});
    }

    MXDict nlp_prob = {
        {"f", cost_fn_},
        {"x", opt_var},
        {"p", opt_dec},
        {"g", g_}
    };

    std::string solver_name = "ipopt";
    Dict nlp_opts;
    nlp_opts["expand"] = true;
    nlp_opts["ipopt.max_iter"] = 5000;
    nlp_opts["ipopt.print_level"] = 0;
    nlp_opts["print_time"] = 0;
    nlp_opts["ipopt.acceptable_tol"] =  1e-6;
    nlp_opts["ipopt.acceptable_obj_change_tol"] = 1e-4;

    solver_ = nlpsol("nlpsol", solver_name, nlp_prob, nlp_opts);
}

void AUTO_OMNI::set_boundary(std::vector<double> x_min, std::vector<double> x_max,
                             std::vector<double> u_min, std::vector<double> u_max)
{
    for (int k = 0; k < prediction_horizons_+1; k++)
    {
        lbx_.push_back(x_min[0]);
        lbx_.push_back(x_min[1]);
        lbx_.push_back(x_min[2]);

        ubx_.push_back(x_max[0]);
        ubx_.push_back(x_max[1]);
        ubx_.push_back(x_max[2]);
    }

    for (int k = 0; k < prediction_horizons_; k++)
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

    // args_ = {
    //     {"lbx", lbx_},
    //     {"ubx", ubx_},
    //     {"lbg", lbg_},
    //     {"ubg", ubg_},
    // };

}

void AUTO_OMNI::input_trajectory(
        Eigen::Vector3d current_states, Eigen::Vector4d current_controls,
        Eigen::Vector3d goal_states, Eigen::Vector4d goal_controls)
{
    // Params

    std::vector<double> variables;
    std::vector<double> params;

    // 3 * (num_states+1) + 4 * num_control

    for (int j = 0; j < prediction_horizons_+1; j++)
    {
         for (int k = 0; k < 3; k++)
         {
            variables.push_back(current_states[k]);
         }
    }

    for (int j = 0; j < prediction_horizons_; j++)
    {
         for (int k = 0; k < 4; k++)
         {
            variables.push_back(current_controls[k]);
         }
    }

    params.push_back(current_states(0));
    params.push_back(current_states(1));
    params.push_back(current_states(2));

    for (int j = 0; j < prediction_horizons_; j++)
    {
         for (int k = 0; k < 3; k++)
         {
            params.push_back(goal_states[k]);
         }
    }

    for (int j = 0; j < prediction_horizons_; j++)
    {
         for (int k = 0; k < 4; k++)
         {
            params.push_back(goal_controls[k]);
         }
    }

    args_ = {
        {"lbx", lbx_},
        {"ubx", ubx_},
        {"lbg", lbg_},
        {"ubg", ubg_},
        {"x0", variables},
        {"p", params}
    };

    // std::cout << args_["x0"] << std::endl;

}

std::vector<double> AUTO_OMNI::optimal_solution()
{
    results_ = solver_(args_);

    // std::cout << "Stage 8" << std::endl;

    std::vector<double> results_all(results_.at("x"));

    // std::cout << result_u[0] << std::endl;

    return results_all;
}

Eigen::VectorXd AUTO_OMNI::forward_kinematic(
    double u1, double u2, double u3, double u4, double theta
)
{
    auto J = Eigen::MatrixXd(3, 4);
    J << sin(a1_+theta), -sin(a2_+theta), sin(a3_+theta), -sin(a4_+theta),
         cos(a1_+theta), -cos(a2_+theta), cos(a3_+theta), -cos(a4_+theta),
         1/(2*L_), 1/(2*L_), 1/(2*L_), 1/(2*L_);

    J = (r_/2)*J;

    auto inv_vec = Eigen::Vector4d(4);

    inv_vec << u1, u2, u3, u4;

    auto for_vec = J*inv_vec;

    return for_vec;
}       

Eigen::VectorXd AUTO_OMNI::inverse_kinematic(
    double vx, double vy, double vyaw, double theta
)
{
    auto J = Eigen::MatrixXd(4, 3);
    J << sin(a1_+theta), cos(a1_+theta), L_,
         sin(a2_+theta), cos(a2_+theta), L_,
         sin(a3_+theta), cos(a3_+theta), L_,
         sin(a4_+theta), cos(a4_+theta), L_;

    J = (1/r_)*J;

    auto for_vec = Eigen::Vector4d(3);

    auto inv_vec = J*for_vec;


    return inv_vec;
}