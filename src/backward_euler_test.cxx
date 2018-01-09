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

TEST(BackwardEulerTest, SimpleLinearIncrease) {
    int digits = 15;

    std::vector<real_t> rates_list = {1, 10, 50};
    std::vector<real_t> dt_list = {0.01, 0.1, 1, 10, 50};
    int nsteps = 2;
    int nentries = 10;
    for (const auto& dt : dt_list) {
        for (const auto& rate : rates_list) {
            colvec_t u(nentries);
            u.fill(0.0);
            auto transform = ConstantShift(rate, nentries);
            auto integrator = BackwardEuler(&transform, u, dt);

            for (int k = 0; k < nsteps; k++) {
                integrator.take_timestep();
                auto state = integrator.get_state();
                for (int l = 0; l < nentries; l++) {
                    EXPECT_NEAR_DIGITS(state[l], (k + 1) * dt * rate, digits);
                }
            }
        }
    }
}
 
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
