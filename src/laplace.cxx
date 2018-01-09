#include <operators/laplace.hxx>

LaplaceOperator1D::LaplaceOperator1D(const int npoints, const real_t right_bc, 
                                     const real_t scale_factor) {
    //compute the discretized difference
    m_dx = 1.0 / (npoints + 1);

    m_laplacian_matrix.resize(npoints, npoints);
    //make the discretized Laplace operator
    for (int i = 0; i < npoints; i++) {
        for (int j = 0; j < npoints; j++) {
            m_laplacian_matrix(i, j) = 0;
        }
    }

    //compute the scale factor but don't multiply it into the laplacian yet
    m_scale_factor = scale_factor * std::pow(m_dx, -2);

    //construct the appropriate discretization of the laplacian matrix
    for (int k = 0; k < npoints; k++) {
        m_laplacian_matrix(k, k) = -2;
        if (k != npoints - 1) {                                  
            m_laplacian_matrix(k, k+1) = 1;
            m_laplacian_matrix(k+1, k) = 1;
        }
    }

    //make the boundary corrector
    m_boundary_term.resize(npoints);
    m_boundary_term(0) = 1.0;
    m_boundary_term(npoints - 1) = right_bc;
    for (int k = 1; k < npoints - 1; k++) {
        m_boundary_term(k) = 0;
    }
}

colvec_t LaplaceOperator1D::apply(const colvec_t& x) const {
    //scale the result at the end to avoid precision loss
    auto unscaled_result = m_laplacian_matrix * x + m_boundary_term;
    return unscaled_result * m_scale_factor;
}

colvec_t LaplaceOperator1D::scale(const colvec_t& x) const {
    return m_laplacian_matrix * m_scale_factor * x;
}

colvec_t LaplaceOperator1D::shift(const colvec_t& x) const {
    return x + m_scale_factor * m_boundary_term;
}

void LaplaceOperator1D::apply(colvec_t& x) const {
    //scale the result at the end to avoid precision loss
    auto unscaled_result = m_laplacian_matrix * x + m_boundary_term;
    auto scaled_result = m_scale_factor * unscaled_result;
    x = scaled_result;
}

void LaplaceOperator1D::scale(colvec_t& x) const {
    x = m_scale_factor * m_laplacian_matrix * x;
}

void LaplaceOperator1D::shift(colvec_t& x) const {
    x += m_boundary_term * m_scale_factor;
}
mat_t LaplaceOperator1D::get_scale() const {
    return m_laplacian_matrix * m_scale_factor;
}

colvec_t LaplaceOperator1D::get_shift() const {
    return m_boundary_term * m_scale_factor;
}

const mat_t& LaplaceOperator1D::get_laplacian() const {
    return m_laplacian_matrix;
}

const colvec_t& LaplaceOperator1D::get_boundary_term() const {
    return m_boundary_term;
}
