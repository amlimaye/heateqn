#include "laplace.hxx"
#include <iostream>

LaplaceOperator1D::LaplaceOperator1D(const int a_npoints, const real_t a_right_bc) {
    //compute the discretized difference
    m_dx = 1.0 / (a_npoints + 1);

    m_laplacian_matrix.resize(a_npoints, a_npoints);
    //make the discretized Laplace operator
    for (int i = 0; i < a_npoints; i++) {
        for (int j = 0; j < a_npoints; j++) {
            m_laplacian_matrix(i, j) = 0;
        }
    }

    real_t scale = std::pow(m_dx, -2);

    for (int k = 0; k < a_npoints; k++) {
        m_laplacian_matrix(k, k) = -2*scale;            //-2 along the diagonal for all indices
        if (k != a_npoints - 1) {                       //if not at the last entry, place 1 to the bottom and right
            m_laplacian_matrix(k, k+1) = 1*scale;
            m_laplacian_matrix(k+1, k) = 1*scale;
        }
    }

    //make the boundary corrector
    m_boundary_term.resize(a_npoints);
    m_boundary_term(0) = scale;
    m_boundary_term(a_npoints - 1) = a_right_bc * scale;
}

colvec_t LaplaceOperator1D::operator()(const colvec_t& a_apply_to) {
    return m_laplacian_matrix * a_apply_to + m_boundary_term;
}
