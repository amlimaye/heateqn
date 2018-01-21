#pragma once
#include <types.hxx>

class Domain2D {
public:
    Domain2D(int npoints_per_dim, real_t value = 0.0);
    
    colvec_t& get_interior_data() {return m_interior_data;};
    colvec_t& get_boundary_data() {return m_boundary_data;};
    int get_npoints_per_dim() const {return m_npoints_per_dim;};

    colvec_t get_boundary_corrector() const;

    void map_over_domain(const std::function<real_t(real_t, real_t)>& f);

private:
    int m_npoints_per_dim;
    colvec_t m_interior_data;
    colvec_t m_boundary_data;
};
