#include <gtest/gtest.h>
#include <math.h>
#include <operators/laplace.hxx>
#include <operators/constant_shift.hxx>
#include <integrators/forward_euler.hxx>

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

/*  computes the stability threshold for the forward Euler integrator, which is
 *  derived to be: dt < [(dx)^2 / (2 * alpha)] sin^2({\pi / 2} * {npts / (npts + 1)})
 *  from the maximum magnitude eigenvalue of the Laplacian
 */
double compute_stability_threshold(real_t alpha, int npoints) {
    real_t dx = 1.0 / (npoints + 1);

    real_t dxsq_over_two_alpha = std::pow(dx, 2) / (2 * alpha);
    real_t n_over_nplusone = (double) npoints / (npoints + 1);
    real_t pi_over_two = M_PI / 2;
    real_t sine_factor = std::pow(std::sin(pi_over_two * n_over_nplusone), 2);

    real_t stability_threshold = dxsq_over_two_alpha / sine_factor;
    return stability_threshold;
}

TEST(ForwardEulerTest, SimpleLinearIncrease) {
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
            auto integrator = ForwardEuler(&transform, u, dt);

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

TEST(ForwardEulerTest, LaplaceStabilityThreshold) {
    std::vector<real_t> alpha_list = {0.5, 1.0, 10.0};
    std::vector<int> npoints_list = {10, 50};
    real_t right_bc = 0.5;

    for (const auto& alpha : alpha_list) {
        for (const auto& npoints : npoints_list) {
            int nsteps = (int) 50000.0 / alpha;
            colvec_t initial_condition;
            make_expo_ic(npoints, right_bc, initial_condition);

            auto laplace = LaplaceOperator1D(npoints, right_bc, alpha);

            auto thresh_dt = compute_stability_threshold(alpha, npoints);

            auto stable_integrator = ForwardEuler(&laplace, initial_condition,
                    0.98 * thresh_dt);
            auto unstable_integrator = ForwardEuler(&laplace, initial_condition,
                    1.02 * thresh_dt);

            for (int t = 0; t < nsteps; t++) {
                stable_integrator.take_timestep();
                unstable_integrator.take_timestep();
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
}
 
TEST(ForwardEulerTest, LongTimeLaplaceConvergence) {
    //TODO: why am I losing digits here?
    int digits = 12;
    std::vector<int> npoints_list = {10, 50, 100};
    std::vector<real_t> alpha_list = {0.1, 1.0, 5.0};

    real_t right_bc = 0.5;

    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            //make the exponential initial condition and the linear steady state vector
            colvec_t initial_condition;
            colvec_t steady_state;
            make_expo_ic(npoints, right_bc, initial_condition);
            make_linear_ic(npoints, right_bc, steady_state);

            //get a dt that leads to stable integration and calculate the # of steps
            real_t dt = 0.99 * compute_stability_threshold(alpha, npoints);
            auto converged_time = 4.0 / alpha;
            int converged_steps = converged_time / dt;

            //integrate for required # of steps
            auto laplace = LaplaceOperator1D(npoints, right_bc, alpha);
            auto integrator = ForwardEuler(&laplace, initial_condition, dt);
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
