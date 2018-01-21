#include <types.hxx>

class Domain2D {
public:
    Domain2D(int npoints_per_dim, real_t value = 0.0);
    
    colvec_t& get_interior_data();
    colvec_t& get_boundary_data();
    colvec_t get_boundary_corrector();

    template<typename Functor>
    void map_over_domain(Functor& f);

private:
    int m_npoints_per_dim;
    colvec_t m_interior_data;
    colvec_t m_boundary_data;
};

#include <domain.txx>
