Domain2D::Domain2D(int npoints_per_dim, real_t value) {
    m_npoints_per_dim = npoints_per_dim;

    m_interior_data.resize(npoints_per_dim * npoints_per_dim);
    m_boundary_data.resize(npoints_per_dim * 4);

    for (int k = 0; k < npoints_per_dim * npoints_per_dim; k++) {
        m_interior_data[k] = value;
    }

    for (int k = 0; k < npoints_per_dim * 4; k++) {
        m_boundary_data[k] = value;
    }
}

colvec_t& Domain2D::get_interior_data() {
    return m_interior_data;
}

colvec_t& Domain2D::get_boundary_data() {
    return m_boundary_data;
}

colvec_t Domain2D::get_boundary_corrector() {
    colvec_t boundary_corrector;
    boundary_corrector.resize(m_npoints_per_dim * m_npoints_per_dim);
    boundary_corrector.setZero();

    //convert x,y index pair into a flat, x-major indexing scheme
    auto get_gidx = [&] (const int x_idx, const int y_idx) -> int {
        return x_idx + (m_npoints_per_dim * y_idx);
    };

    for (int k = 0; k < m_npoints_per_dim; k++) {
        //left
        boundary_corrector(get_gidx(0, k)) += m_boundary_data(k);

        //bottom
        boundary_corrector(get_gidx(k, 0)) += 
            m_boundary_data(m_npoints_per_dim + k);

        //top
        boundary_corrector(get_gidx(m_npoints_per_dim - 1, k)) += 
            m_boundary_data(m_npoints_per_dim * 2 + k);

        //right
        boundary_corrector(get_gidx(k, m_npoints_per_dim - 1)) += 
            m_boundary_data(m_npoints_per_dim * 3 + k);
    }

    return boundary_corrector;
}

template <typename Functor>
void Domain2D::map_over_domain(Functor& f) {
    //convert x,y index pair into a flat, x-major indexing scheme
    auto get_gidx = [&] (const int x_idx, const int y_idx) -> int {
        return x_idx + (m_npoints_per_dim * y_idx);
    };

    for (int k = 0; k < m_npoints_per_dim; k++) {
        double x = (double) (k+1) / (m_npoints_per_dim+1);
        for (int l = 0; l < m_npoints_per_dim; l++) {
            double y = (double) (l+1) / (m_npoints_per_dim+1);

            m_interior_data(get_gidx(k, l)) = f(x, y);
        }
    }

    //left + right boundaries
    for (int k = 0; k < m_npoints_per_dim; k++) {
        double y = (double) (k+1) / (m_npoints_per_dim+1);
        m_boundary_data(k) = f(0.0, y);
        m_boundary_data(m_npoints_per_dim * 2 + k) = f(1.0, y);
    }

    //top + bottom boundaries
    for (int k = 0; k < m_npoints_per_dim; k++) {
        double x = (double) (k+1) / (m_npoints_per_dim+1);
        m_boundary_data(m_npoints_per_dim + k) =  f(x, 0.0);
        m_boundary_data(m_npoints_per_dim * 3 + k) =  f(x, 1.0);
    }
}
