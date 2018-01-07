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

    //compute the scale factor but don't multiply it into the laplacian yet
    m_scale = std::pow(m_dx, -2);

    //construct the appropriate discretization of the laplacian matrix
    for (int k = 0; k < a_npoints; k++) {
        m_laplacian_matrix(k, k) = -2;
        if (k != a_npoints - 1) {                                  
            m_laplacian_matrix(k, k+1) = 1;
            m_laplacian_matrix(k+1, k) = 1;
        }
    }

    //make the boundary corrector
    m_boundary_term.resize(a_npoints);
    m_boundary_term(0) = 1.0;
    m_boundary_term(a_npoints - 1) = a_right_bc;
    for (int k = 1; k < a_npoints - 1; k++) {
        m_boundary_term(k) = 0;
    }
}

colvec_t LaplaceOperator1D::apply(const colvec_t& a_apply_to) {
    //scale the result before returning to avoid precision loss
    auto unscaled_result = m_laplacian_matrix * a_apply_to + m_boundary_term;
    return unscaled_result * m_scale;
}

const mat_t& LaplaceOperator1D::get_laplacian() const {
    return m_laplacian_matrix;
}

const mat_t LaplaceOperator1D::get_scaled_laplacian() const {
    return m_scale*m_laplacian_matrix;
}

const colvec_t& LaplaceOperator1D::get_boundary_term() const {
    return m_boundary_term;
}

const colvec_t LaplaceOperator1D::get_scaled_boundary_term() const {
    return m_scale*m_boundary_term;
}
