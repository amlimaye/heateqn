#include "types.hxx"

class LaplaceOperator1D {
public:
    //ctor
    LaplaceOperator1D(const int a_npoints, const real_t a_right_bc);

    //apply to an column vector
    colvec_t apply(const colvec_t& a_apply_to);

    //accessors
    real_t              get_dx() const {return m_dx;};
    const mat_t&        get_laplacian() const;
    const colvec_t&     get_boundary_term() const;

    //Warning: these make copies!
    const mat_t         get_scaled_laplacian() const;
    const colvec_t      get_scaled_boundary_term() const;

private:
    real_t      m_dx;
    real_t      m_scale;
    mat_t       m_laplacian_matrix;
    colvec_t    m_boundary_term;
};
