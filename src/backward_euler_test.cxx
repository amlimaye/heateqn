#include <gtest/gtest.h>
#include <math.h>
#include "laplace.hxx"
#include "constant_shift.hxx"
#include "backward_euler.hxx"

#define EXPECT_NEAR_DIGITS(x,y,d) EXPECT_NEAR(x, y, std::pow(10, -(d)))

void make_linear_ic(int npoints, const real_t& right_bc, colvec_t& u) {
    u.resize(npoints);
    auto drop_per_interval = (1.0 - right_bc) / (npoints + 1);
    for (int l = 0; l < npoints; l++) {
        u[l] = 1.0 - (l+1)*drop_per_interval;
    }
};

void make_expo_ic(int npoints, const real_t& right_bc, colvec_t& u) {
    u.resize(npoints);
    double decay_rate = -1.0 * std::log(right_bc);

    for (int k = 0; k < npoints; k++) {
        double x = (double) (k+1) / (npoints+1);
        u(k) = std::exp(-1.0 * decay_rate * x);
    } 
};


/*TEST(ForwardEulerTest, LaplaceStabilityThreshold) {
    std::vector<real_t> alpha_list = {0.5, 1.0, 10.0};
    std::vector<int> npoints_list = {10, 50, 100};
    real_t right_bc = 0.5;

    for (const auto& alpha : alpha_list) {
        for (const auto& npoints : npoints_list) {
            int nsteps = (int) 50000.0 / alpha;
            colvec_t initial_condition;
            make_expo_ic(npoints, right_bc, initial_condition);

            auto laplace = LaplaceOperator1D(npoints, right_bc);

            auto thresh_dt = compute_stability_threshold(alpha, npoints);
            auto stable_integrator = ForwardEuler(initial_condition, 0.98 *
                    thresh_dt);
            auto unstable_integrator = ForwardEuler(initial_condition, 1.02 *
                    thresh_dt);

            for (int t = 0; t < nsteps; t++) {
                stable_integrator.take_timestep(alpha *
                        laplace.apply(stable_integrator.get_state()));
                unstable_integrator.take_timestep(alpha *
                        laplace.apply(unstable_integrator.get_state()));
            }

            //stable integration should always produce values of magnitude less
            //than 1.0 because we scaled the initial condition to a max value
            //of 1.0.  unstable integration, on the other hand, will eventually
            //produce values that get very large (they may just NaN out)
            auto stable_state = stable_integrator.get_state();
            auto unstable_state = unstable_integrator.get_state();
            for (int k = 0; k < npoints; k++) {
                EXPECT_LE(std::abs(stable_state[k]), 1.0);
                EXPECT_TRUE(std::abs(unstable_state[k]) > 1.0 ||
                        isnan(unstable_state[k]));
            }
        }
    }
}*/
 
TEST(BackwardEulerTest, LongTimeLaplaceConvergence) {
    //TODO: why am I losing digits here?
    int digits = 12;
    std::vector<int> npoints_list = {10};
    std::vector<real_t> alpha_list = {1.0};

    real_t right_bc = 0.5;

    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            //make the exponential initial condition and the linear steady state vector
            colvec_t initial_condition;
            colvec_t steady_state;
            make_expo_ic(npoints, right_bc, initial_condition);
            make_linear_ic(npoints, right_bc, steady_state);

            //get a dt that leads to stable integration and calculate the # of steps
            real_t dt = 5.00;
            auto converged_time = 1000.0 / alpha;
            int converged_steps = converged_time / dt;

            //integrate for required # of steps
            auto laplace = LaplaceOperator1D(npoints, right_bc);
            auto integrator = BackwardEuler(&laplace, initial_condition, dt);
            for (int t = 0; t < converged_steps; t++) {
                integrator.take_timestep();
            }
            auto state = integrator.get_state();

            //check if we are near the linear steady state
            for (int k = 0; k < npoints; k++) {
                EXPECT_NEAR_DIGITS(state[k], steady_state[k], digits);
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
