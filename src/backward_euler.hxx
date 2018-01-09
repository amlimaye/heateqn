#include "types.hxx"
#include "affine_transform.hxx"

class BackwardEuler {
public:
    BackwardEuler (AffineTransform* a_transform, colvec_t a_ic, real_t a_dt);

    void                take_timestep();
    const colvec_t&     get_state() const;
    const real_t&       get_dt() const;

protected:
    colvec_t            m_curr_state;
    real_t              m_dt;
    Eigen::LLT<mat_t>   m_llt_decomposed;
    colvec_t            m_rhs_shift;
};
