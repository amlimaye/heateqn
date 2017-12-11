#include "types.hxx"

class LaplaceOperator1D {
public:
    LaplaceOperator1D(const int a_npoints, const real_t a_right_bc);
    colvec_t operator()(const colvec_t& a_apply_to);

private:
    real_t      m_dx;
    mat_t       m_laplacian_matrix;
    colvec_t    m_boundary_term;
};
