#include "backward_euler.hxx"

BackwardEuler::BackwardEuler(AffineTransform* a_transform, colvec_t a_ic, 
                             real_t a_dt) {
    m_curr_state = a_ic;
    m_dt = a_dt;

    //cache an LLT decomposition of (I - dt * L) for doing matrix solves
    //LLT decomposition requires positive definite-ness, which works for the
    //laplacian, but this will need to be changed later for generalization
    //purposes
    mat_t identity;
    mat_t A;

    auto scale_matrix = a_transform->get_scale();

    identity.setIdentity(scale_matrix.rows(), scale_matrix.cols());
    A = identity - m_dt * scale_matrix;
    m_llt_decomposed = A.llt();

    m_rhs_shift = m_dt * a_transform->get_shift();
};

void BackwardEuler::take_timestep() {
    m_curr_state = m_llt_decomposed.solve(m_curr_state + m_rhs_shift);
}

const colvec_t& BackwardEuler::get_state() const {
    return m_curr_state;
}

const real_t& BackwardEuler::get_dt() const {
    return m_dt;
}
