#include "forward_euler.hxx"

ForwardEuler::ForwardEuler(colvec_t a_ic, real_t a_dt) :
    m_curr_state(a_ic),
    m_dt(a_dt)
{};

void ForwardEuler::take_timestep(const colvec_t& a_rhs) {
    m_curr_state += a_rhs * m_dt;
}

const colvec_t& ForwardEuler::get_state() const {
    return m_curr_state;
}

const real_t& ForwardEuler::get_dt() const {
    return m_dt;
}
