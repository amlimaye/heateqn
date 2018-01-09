#include "types.hxx"
#include "affine_transform.hxx"

class ForwardEuler {
public:
    ForwardEuler(AffineTransform* a_transform, colvec_t a_ic, real_t a_dt);

    void                take_timestep();
    const colvec_t&     get_state() const;
    const real_t&       get_dt() const;

protected:
    AffineTransform*    m_transform;
    colvec_t            m_curr_state;
    real_t              m_dt;
};
