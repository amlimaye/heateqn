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

int main() {
    constexpr real_t dt = 0.005;
    constexpr real_t right_bc = 0.15;
    constexpr real_t alpha = 0.50;
    constexpr int npoints = 10;
    constexpr int nsteps = 60;

    colvec_t ic = make_expo_ic(npoints, right_bc);
    Timestepper ts(ic, dt);
    LaplaceOperator1D laplace(npoints, right_bc);

    std::cout << ts.get_state() << std::endl;
    for (int t = 0; t < nsteps; t++) {
        ts.take_timestep(alpha * laplace(ts.get_state()));

        if ((t+1) % 10 == 0) {
            std::cout << ts.get_state() << std::endl;
        }
    }

    return 0;
}
