#ifndef AUTO_MECANUM_HPP__
#define AUTO_MECANUM_HPP__

#include <iostream>
#include <Eigen/Dense>
#include <casadi/casadi.hpp>
#include <map>
#include <chrono>

using namespace casadi;

class AUTO_MECANUM
{
    private:
        // Mecanum and MPC params
        double Lx_;
        double Ly_;
        double wheel_radius_;
        double dt_;
        double sim_time_;
        int prediction_horizon_;

        Slice all;
        Slice slice_states;
        Slice slice_controls;
        DM Q_ = DM::zeros(3,3);
        DM R_ = DM::zeros(4,4);
        MX rot_mat = MX::ones(3,3);
        MX J_for = MX::ones(3,4);
        MX cost_fn = 0.0;

        MX g;

        std::vector<double> lbx_;
        std::vector<double> ubx_;

        DM lbg_ = 0.0;
        DM ubg_ = 0.0;

        DMDict args_;
        Function solver_;

        // std::vector<MX> RHS;

        // std::vector<MX> RHS;
        MX RHS;

        // MAP
        std::map<std::string, DM> results_;

        // Casadi Function
        Function system_kinematic_;

        // Eigen States
        Eigen::Vector3d current_states;
        Eigen::Vector4d current_controls;

        template <typename T1, typename T2>
        T1 RK4_function(T1 state, T1 control, T2 dynamics)
        {
            std::vector<T1> input;
            input.push_back(state);
            input.push_back(control);
            T1 k1 = dynamics(input).at(0);
            input[0] = state + dt_ / 2.0 * k1;
            T1 k2 = dynamics(input).at(0);
            input[0] = state + dt_ / 2.0 * k2;
            T1 k3 = dynamics(input).at(0);
            input[0] = state + dt_ * k3;
            T1 k4 = dynamics(input).at(0);
            return state + dt_ / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }

    public:

        AUTO_MECANUM(
            double Lx, double Ly, double wheel_radius,
            double dt, double sim_time, int prediction_horizon
        );

        virtual ~AUTO_MECANUM();

        Eigen::VectorXd forward_kinematic(double theta, double w1, double w2, double w3, double w4);
        Eigen::VectorXd inverse_kinematic(double theta, double vx, double vy, double vth);
        Eigen::VectorXd shift_timestep(double dt, double &t, std::vector<double> x0, std::vector<double> &u);

        void setup_auto_mecanum();

        void set_boundary(std::vector<double> x_min, std::vector<double> x_max,
                          std::vector<double> u_min, std::vector<double> u_max);

        void input_trajectory(Eigen::Vector3d current_states, Eigen::Vector4d current_controls,
                              Eigen::Vector3d goal_states, Eigen::Vector4d goal_controls);

        std::vector<double> optimal_solution();


};

#endif