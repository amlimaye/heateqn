#include <vector>
#include <tuple>
#include <iostream>
#include <integrators/forward_euler.hxx>
#include <integrators/backward_euler.hxx>
#include <operators/scale.hxx>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

real_t exact_solution(const real_t initial, const real_t alpha, 
                      const real_t time) {
    return initial * std::exp(-alpha * time);
}

typedef std::vector<std::tuple<real_t, real_t>> dts_and_gerrs_t;
typedef std::vector<std::tuple<real_t, dts_and_gerrs_t>> to_python_t;

to_python_t get_global_errors_forward_euler() {
    //not dealing with a vector-valued y here. just one entry is enough
    int npoints = 1;
    real_t initial = 1.0;
    colvec_t ic(npoints);
    for (int k = 0; k < npoints; k++) {
        ic(k) = initial;
    }

    //set the rate and compute the final time to be 4 time constants
    real_t rate = 2.0;

    auto inv_rate = 1.0 / rate;
    std::vector<real_t> final_times_list = {inv_rate, 4.0 * inv_rate, 
                                            6.0 * inv_rate, 8.0 * inv_rate};

    //set up the RHS to be the scaling operator, results in expo. decay
    auto op = Scale(npoints, -rate);

    //make list of dt values
    std::vector<real_t> dt_list = {0.001, 0.002, 0.005, 0.01, 0.02, 0.04, 0.05,
                                   0.1, 0.2};

    //compute the global error for each dt
    to_python_t out;
    for (const auto& final_time : final_times_list) {
        dts_and_gerrs_t dts_and_gerrs;
        for(const auto& dt : dt_list) {
            int nsteps = (int) std::ceil(final_time / dt);
            auto integrator = ForwardEuler(&op, ic, dt);

            real_t time = 0.0; 
            for (int t = 0; t < nsteps; t++) {
                integrator.take_timestep();
                time += dt;
            }

            auto state = integrator.get_state();
            auto global_error = std::abs(state(0) - 
                                         exact_solution(initial, rate, time));

            dts_and_gerrs.push_back(std::make_tuple(dt, global_error));
        }
        out.push_back(std::make_tuple(final_time, dts_and_gerrs));
    }

    return out;
}

to_python_t get_global_errors_backward_euler() {
    //not dealing with a vector-valued y here. just one entry is enough
    int npoints = 1;
    real_t initial = 1.0;
    colvec_t ic(npoints);
    for (int k = 0; k < npoints; k++) {
        ic(k) = initial;
    }

    //set the rate and compute the final time to be 4 time constants
    real_t rate = 2.0;

    auto inv_rate = 1.0 / rate;
    std::vector<real_t> final_times_list = {inv_rate, 4.0 * inv_rate, 
                                            6.0 * inv_rate, 8.0 * inv_rate};

    //set up the RHS to be the scaling operator, results in expo. decay
    auto op = Scale(npoints, -rate);

    //make list of dt values
    std::vector<real_t> dt_list = {0.001, 0.002, 0.005, 0.01, 0.02, 0.04, 0.05,
                                   0.1, 0.2, 0.50, 1.0, 2.0};

    //compute the global error for each dt
    to_python_t out;
    for (const auto& final_time : final_times_list) {
        dts_and_gerrs_t dts_and_gerrs;
        for(const auto& dt : dt_list) {
            int nsteps = (int) std::ceil(final_time / dt);
            auto integrator = BackwardEuler(&op, ic, dt);

            real_t time = 0.0; 
            for (int t = 0; t < nsteps; t++) {
                integrator.take_timestep();
                time += dt;
            }

            auto state = integrator.get_state();
            auto global_error = std::abs(state(0) - 
                                         exact_solution(initial, rate, time));

            dts_and_gerrs.push_back(std::make_tuple(dt, global_error));
        }
        out.push_back(std::make_tuple(final_time, dts_and_gerrs));
    }

    return out;
}


namespace py = pybind11;

PYBIND11_MODULE(pyexample_convergence_rate, m) {
    m.def("forward_euler", &get_global_errors_forward_euler);
    m.def("backward_euler", &get_global_errors_backward_euler);
}
