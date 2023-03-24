#ifndef __MPC_MECANUM_HPP__
#define __MPC_MECANUM_HPP__
#include <Eigen/Dense>
#include <chrono>
#include <casadi/casadi.hpp>
#include <iostream>

namespace ca = casadi;

class MPCMECANUM
{
    public:
        MPCMECANUM(
            const double Lx_, const double Ly_, const double wheel_radius_,
            const double time_step_, const double sim_time_, int prediction_horizons_
        );

        virtual ~MPCMECANUM();

        void setup_mpc();
        void set_boundary(Eigen::VectorXd u_min, Eigen::VectorXd u_max,
                          Eigen::VectorXd x_min, Eigen::VectorXd x_max);

        void trajectory_generation(Eigen::VectorXd init_states, Eigen::VectorXd target_states,
                                   Eigen::VectorXd init_controls, Eigen::VectorXd target_controls);

        std::vector<double> get_optimal_solution();

        // Function get_dynamic_model(MX states, MX controls);
        // ---- Eigen Kinematic Model ---- //
        Eigen::VectorXd kinematic_eigen(Eigen::VectorXd states, Eigen::VectorXd controls);

        double forward_kinematic();
        double inverse_kinematic();


    private:
        // ---- Mecanum Parameters ---- //
        double Lx_;
        double Ly_;
        double wheel_radius_;

        // ---- MPC Parameters ---- //
        int prediction_horizon_;
        double time_step_;
        double sim_time_;
        int num_states_ = 3;
        int num_controls_ = 4;
        // ---- MPC Constraint and Cost Function
        ca::MX cost_fn_;
        ca::MX g_;

        std::vector<double> x_min_;
        std::vector<double> x_max_;
        std::vector<double> u_min_;
        std::vector<double> u_max_;

        // ---- MPC Argument ---- //
        ca::DMDict args_;
        // ---- MPC Index ---- //
        ca::Slice all;
        std::vector<ca::MX> g_;


        // ---- Ipopt solver ---- //
        ca::Function system_kinematic_model();

        ca::Function solver_;

        std::map<std::string, ca::DM> opt_result_;

        template<typename T1, typename T2>
        T1 get_kinematic_model(T1 states, T2 controls)
        {
            ca::DM rot_mat = DM::zeros(3,3);
            rot_mat(0,0) = ca::cos(states[2]);
            rot_mat(0,1) = ca::sin(states[2]);
            rot_mat(0,2) = 0;
            rot_mat(1,0) = -ca::sin(states[2]);
            rot_mat(1,1) = ca::cos(states[2]);
            rot_mat(1,2) = 0;
            rot_mat(2,0) = 0;
            rot_mat(2,1) = 0;
            rot_mat(2,2) = 1;

            ca::DM J_for = DM::zeros(3,4);
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

            ca::DM rhs(3,1);

            rhs = ca::MX::mtimes({rot_mat.T(), J_for, controls});

            return rhs;
        }

};


#endif