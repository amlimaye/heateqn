#include "backward_euler.hxx"

BackwardEuler::BackwardEuler(colvec_t a_ic, real_t a_dt) :
    m_curr_state(a_ic),
    m_dt(a_dt)
{};

void BackwardEuler::take_timestep(const mat_t& a_operator, const colvec_t& a_boundary_term) {
    mat_t identity;
    identity.setIdentity(a_operator.rows(), a_operator.cols());

    mat_t A = identity - m_dt * a_operator;

    //LLT decomposition requires positive definite-ness, which works for the
    //laplacian, but this will need to be changed later for generalization
    //purposes
    m_curr_state = A.llt().solve(m_curr_state + m_dt*a_boundary_term);
}

const colvec_t& BackwardEuler::get_state() const {
    return m_curr_state;
}

const real_t& BackwardEuler::get_dt() const {
    return m_dt;
}
