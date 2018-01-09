#include "types.hxx"

class BackwardEuler {
public:
    BackwardEuler (colvec_t a_ic, real_t a_dt);

    void                take_timestep(const mat_t& a_operator, const colvec_t& a_boundary_term);
    const colvec_t&     get_state() const;
    const real_t&       get_dt() const;

protected:
    colvec_t m_curr_state;
    real_t m_dt;
};
