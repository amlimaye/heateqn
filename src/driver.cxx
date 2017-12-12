#include "timestepper.hxx"
#include "laplace.hxx"
#include <iostream>

colvec_t make_expo_ic(const int npoints, const real_t right_bc) {
    //initialize the output vector
    colvec_t out(npoints);

    //compute the decay rate of the exponential
    double decay_rate = -1.0 * std::log(right_bc);

    //assign values sampled from this exponential to the initial condition
    for (int k = 0; k < npoints; k++) {
        double x = (double) (k+1) / (npoints+1);
        out(k) = std::exp(-1.0 * decay_rate * x);
    }

    return out;
}

/*  computes the stability threshold for the forward Euler integrator, which is
 *  derived to be: dt < [(dx)^2 / (2 * alpha)] sin^2({\pi / 2} * {npts / (npts + 1)})
 */
double compute_stability_threshold(real_t alpha, int npoints) {
    real_t dx = 1.0 / npoints;

    double dxsq_over_two_alpha = std::pow(dx, 2) / (2 * alpha);
    double n_over_nplusone = (double) npoints / (npoints + 1);
    double pi_over_two = M_PI / 2;
    double sine_factor = std::pow(std::sin(pi_over_two * n_over_nplusone), 2);

    double stability_threshold = dxsq_over_two_alpha / sine_factor;
    return stability_threshold;
}

void fwd_euler_stability_test(const real_t alpha, const int npoints) {
    //set the right boundary condition (arbitrarily) to zero, make the IC
    constexpr real_t right_bc = 0.0;
    colvec_t ic = make_expo_ic(npoints, right_bc);

    //compute the stability threshold
    auto thresh_dt = compute_stability_threshold(alpha, npoints);

    //define the nsteps and dt's for the unstable and stable integrators
    auto thresh_dt_unstable = thresh_dt * 1.1;
    int nsteps_unstable = 4 * alpha / thresh_dt_unstable;

    auto thresh_dt_stable = thresh_dt * 0.8;
    int nsteps_stable = 4 * alpha / thresh_dt_stable;

    //construct the timesteppers
    Timestepper this_will_fail(ic, thresh_dt_unstable);
    Timestepper this_will_work(ic, thresh_dt_stable);

    //construct the laplace operator
    LaplaceOperator1D laplace(npoints, right_bc);

    //actually run the test!
    std::cout << "Stability threshold dt = " << thresh_dt << std::endl;

    std::cout << "Running unstable integrator, dt = " << this_will_fail.get_dt() << std::endl;
    for (int t = 0; t < nsteps_unstable; t++) {
       this_will_fail.take_timestep(alpha * laplace.apply(this_will_fail.get_state()));
    }
    std::cout << "|u|_inf at t_final = (4 * alpha) = " << 4 * alpha << std::endl;
    std::cout << this_will_fail.get_state().maxCoeff() << std::endl;

    std::cout << "Running stable integrator, dt = " << this_will_work.get_dt() << std::endl;
    for (int t = 0; t < nsteps_stable; t++) {
       this_will_work.take_timestep(alpha * laplace.apply(this_will_work.get_state()));
    }
    std::cout << "|u|_inf at t_final = (4 * alpha) = " << 4 * alpha << std::endl;
    std::cout << this_will_work.get_state().maxCoeff() << std::endl;
}

int main() {
    constexpr real_t alpha = 2.00;
    constexpr int npoints = 10;
    fwd_euler_stability_test(alpha, npoints);
}
