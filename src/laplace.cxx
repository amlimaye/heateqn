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

LaplaceOperator2D::LaplaceOperator2D(const int npoints_per_dim, 
                                     const colvec_t& boundary_values,
                                     const real_t scale_factor) {
    //compute the discretized differences
    m_dx = 1.0 / (npoints_per_dim + 1);

    //compute a few convenience variables
    auto npoints = npoints_per_dim * npoints_per_dim;

    //convert x,y index pair into a flat, x-major indexing scheme
    auto get_gidx = [&] (const int x_idx, const int y_idx) -> int {
        return x_idx + (npoints_per_dim * y_idx);
    };

    //zero out the laplacian matrix
    m_laplacian_matrix.setZero(npoints, npoints);

    /* TODO: this could be refactored to take advantage of Eigen's
     * symmetrization view functionality, that way we would have to set less
     * elements... but perhaps this maps to a sparse matrix refactor in the near
     * future a little more elegantly.
     */
    
    //construct the laplacian matrix
    m_laplacian_matrix.resize(npoints, npoints);
    for (int x_idx = 0; x_idx < npoints_per_dim; x_idx++) {
        for (int y_idx = 0; y_idx < npoints_per_dim; y_idx++) {
            //set the diagonal term
            auto diag_idx = get_gidx(x_idx, y_idx);
            m_laplacian_matrix(diag_idx, diag_idx) = -4;

            //set the term to the right
            if (x_idx < npoints_per_dim-1) {
                auto right_idx = get_gidx(x_idx + 1, y_idx);
                m_laplacian_matrix(diag_idx, right_idx) = 1;
            }

            //set the term below
            if (y_idx < npoints_per_dim-1) {
                auto down_idx = get_gidx(x_idx, y_idx + 1);
                m_laplacian_matrix(diag_idx, down_idx) = 1;
            }

            //set the term above
            if (y_idx > 0) {
                auto above_idx = get_gidx(x_idx, y_idx - 1);
                m_laplacian_matrix(diag_idx, above_idx) = 1;
            }

            //set the term to the left
            if (x_idx > 0) {
                auto left_idx = get_gidx(x_idx - 1, y_idx);
                m_laplacian_matrix(diag_idx, left_idx) = 1;
            }
        }
    }

    //compute the scale factor but don't multiply it into the laplacian yet
    m_scale_factor = scale_factor * std::pow(m_dx, -2.0);

    //make the boundary corrector
    m_boundary_term.resize(npoints);
    m_boundary_term.setZero(npoints);
    
    /* boundary_values is a flat array of four different boundary condition
     * vectors, the indexing scheme is illustrated by this diagram:
     *
     *         (4)
     *      ----------
     *     |          |
     *     |          |
     * (1) |          | (3)
     *     |          |
     *     |          |
     *      ----------
     *         (2)
     *
     * each edge of the boundary array must have npoints_per_dim elements, so 
     * in total we must have 4 * npoints_per_dim elements
     */
    assert(boundary_values.size() == 4 * npoints_per_dim);
    
    //left boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        m_boundary_term(get_gidx(0, k)) += 
            boundary_values(k);
    }

    //bottom boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        m_boundary_term(get_gidx(k, 0)) += 
            boundary_values(npoints_per_dim + k);
    }

    //right boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        m_boundary_term(get_gidx(npoints_per_dim - 1, k)) += 
            boundary_values(npoints_per_dim * 2 + k);
    }

    //top boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        m_boundary_term(get_gidx(k, npoints_per_dim - 1)) += 
            boundary_values(npoints_per_dim * 3 + k);
    }
}

colvec_t LaplaceOperator1D::apply(const colvec_t& x) const {
    //scale the result at the end to avoid precision loss
    auto unscaled_result = m_laplacian_matrix * x + m_boundary_term;
    return unscaled_result * m_scale_factor;
}

colvec_t LaplaceOperator2D::apply(const colvec_t& x) const {
    //scale the result at the end to avoid precision loss
    auto unscaled_result = m_laplacian_matrix * x + m_boundary_term;
    return unscaled_result * m_scale_factor;
}

colvec_t LaplaceOperator1D::scale(const colvec_t& x) const {
    return m_laplacian_matrix * m_scale_factor * x;
}

colvec_t LaplaceOperator2D::scale(const colvec_t& x) const {
    return m_laplacian_matrix * m_scale_factor * x;
}

colvec_t LaplaceOperator1D::shift(const colvec_t& x) const {
    return x + m_scale_factor * m_boundary_term;
}

colvec_t LaplaceOperator2D::shift(const colvec_t& x) const {
    return x + m_scale_factor * m_boundary_term;
}

void LaplaceOperator1D::apply(colvec_t& x) const {
    //scale the result at the end to avoid precision loss
    auto unscaled_result = m_laplacian_matrix * x + m_boundary_term;
    auto scaled_result = m_scale_factor * unscaled_result;
    x = scaled_result;
}

void LaplaceOperator2D::apply(colvec_t& x) const {
    //scale the result at the end to avoid precision loss
    auto unscaled_result = m_laplacian_matrix * x + m_boundary_term;
    auto scaled_result = m_scale_factor * unscaled_result;
    x = scaled_result;
}

void LaplaceOperator1D::scale(colvec_t& x) const {
    x = m_scale_factor * m_laplacian_matrix * x;
}

void LaplaceOperator2D::scale(colvec_t& x) const {
    x = m_scale_factor * m_laplacian_matrix * x;
}

void LaplaceOperator1D::shift(colvec_t& x) const {
    x += m_boundary_term * m_scale_factor;
}

void LaplaceOperator2D::shift(colvec_t& x) const {
    x += m_boundary_term * m_scale_factor;
}

mat_t LaplaceOperator1D::get_scale() const {
    return m_laplacian_matrix * m_scale_factor;
}

mat_t LaplaceOperator2D::get_scale() const {
    return m_laplacian_matrix * m_scale_factor;
}

colvec_t LaplaceOperator1D::get_shift() const {
    return m_boundary_term * m_scale_factor;
}

colvec_t LaplaceOperator2D::get_shift() const {
    return m_boundary_term * m_scale_factor;
}

const mat_t& LaplaceOperator1D::get_laplacian() const {
    return m_laplacian_matrix;
}

const mat_t& LaplaceOperator2D::get_laplacian() const {
    return m_laplacian_matrix;
}

const colvec_t& LaplaceOperator1D::get_boundary_term() const {
    return m_boundary_term;
}

const colvec_t& LaplaceOperator2D::get_boundary_term() const {
    return m_boundary_term;
}
