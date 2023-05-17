#ifndef AUTO_OMNI_HPP__
#define AUTO_OMNI_HPP__

#include <memory>
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <Eigen/Dense>
#include <math.h>
#include <casadi/casadi.hpp>

using namespace casadi;
using namespace std::literals;
using std::placeholders::_1;


class AUTO_OMNI
{
    public:
        AUTO_OMNI(
            double L, double r, double dt,
            int prediction_horizons, int n_states, int n_controls
        );

        virtual ~AUTO_OMNI();

        Eigen::VectorXd forward_kinematic(double u1, double u2, double u3, double u4, double theta);
        Eigen::VectorXd inverse_kinematic(double vx, double vy, double vth, double theta);
        

        void setup_mpc();

        void set_boundary(std::vector<double> x_min, std::vector<double> x_max,
                          std::vector<double> u_min, std::vector<double> u_max);

        void input_trajectory(Eigen::Vector3d current_states, Eigen::Vector4d current_controls,
                              Eigen::Vector3d goal_states, Eigen::Vector4d goal_controls);
        
        std::vector<double> optimal_solution();

    private:
        // Omni and MPC params
        double r_, L_, a1_, a2_, a3_, a4_, dt_;

        int prediction_horizons_;
        int n_states_;
        int n_controls_;
        
        Slice all;
        Slice slice_states;
        Slice slice_controls;
        DM Q_ = DM::zeros(3, 3);
        DM R_ = DM::zeros(4, 4);
        MX J_for = MX::ones(3, 4);
        MX cost_fn_ = 0.0;
        MX rhs_;
        MX g_;

        std::vector<double> lbx_;
        std::vector<double> ubx_;

        std::map<std::string, DM> results_;

        DM lbg_;
        DM ubg_;

        DMDict args_;
        Function solver_;
        Function system_kinematic_;

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

};

#endif