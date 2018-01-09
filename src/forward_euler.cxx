#include <integrators/forward_euler.hxx>

ForwardEuler::ForwardEuler(AffineTransform* a_transform, colvec_t a_ic, 
                           real_t a_dt) :
    m_transform(a_transform),
    m_curr_state(a_ic),
    m_dt(a_dt)
{};

void ForwardEuler::take_timestep() {
    m_curr_state += m_dt * 
        m_transform->apply(static_cast<const colvec_t>(m_curr_state));
}

const colvec_t& ForwardEuler::get_state() const {
    return m_curr_state;
}

const real_t& ForwardEuler::get_dt() const {
    return m_dt;
}
