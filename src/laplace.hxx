#include "types.hxx"

class LaplaceOperator1D {
public:
    LaplaceOperator1D(const int a_npoints, const real_t a_right_bc);
    colvec_t apply(const colvec_t& a_apply_to);

    //accessors
    real_t              get_dx() const {return m_dx;};
    const mat_t&        get_laplacian() const {return m_laplacian_matrix;};
    const colvec_t&     get_boundary_term() const {return m_boundary_term;};

private:
    real_t      m_dx;
    mat_t       m_laplacian_matrix;
    colvec_t    m_boundary_term;
};
