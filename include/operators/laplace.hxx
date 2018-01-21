#include <types.hxx>
#include <vector>
#include <operators/affine_transform.hxx>
#include <domain.hxx>

class LaplaceOperator1D : public AffineTransform {
public:
    //ctor
    LaplaceOperator1D(const int npoints, const real_t right_bc,
                      const real_t scale_factor = 1.0);

    colvec_t            apply(const colvec_t& x) const;
    void                apply(colvec_t& x) const;

    colvec_t            scale(const colvec_t& x) const;
    colvec_t            shift(const colvec_t& x) const;

    void                scale(colvec_t& x) const;
    void                shift(colvec_t& x) const;

    mat_t               get_scale() const;
    colvec_t            get_shift() const;

    real_t              get_dx() const {return m_dx;};
    const mat_t&        get_laplacian() const;
    const colvec_t&     get_boundary_term() const;

private:
    real_t      m_dx;
    real_t      m_scale_factor;
    mat_t       m_laplacian_matrix;
    colvec_t    m_boundary_term;
};

class LaplaceOperator2D : public AffineTransform {
public:
    //ctor
    LaplaceOperator2D(const Domain2D& domain, 
                      const real_t scale_factor = 1.0);

    colvec_t            apply(const colvec_t& x) const;
    void                apply(colvec_t& x) const;

    colvec_t            scale(const colvec_t& x) const;
    colvec_t            shift(const colvec_t& x) const;

    void                scale(colvec_t& x) const;
    void                shift(colvec_t& x) const;

    mat_t               get_scale() const;
    colvec_t            get_shift() const;

    real_t              get_dx() const {return m_dx;};
    const mat_t&        get_laplacian() const;
    const colvec_t&     get_boundary_term() const;

private:
    real_t      m_dx;
    real_t      m_scale_factor;
    mat_t       m_laplacian_matrix;
    colvec_t    m_boundary_term;
};
