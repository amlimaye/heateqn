#include <types.hxx>
#include <operators/affine_transform.hxx>

class ConstantShift : public AffineTransform {
public:
    ConstantShift(const int a_nrows, const real_t a_shift);

    colvec_t        apply(const colvec_t& x) const;
    void            apply(colvec_t& x) const;

    colvec_t        scale(const colvec_t& x) const;
    colvec_t        shift(const colvec_t& x) const;

    void            scale(colvec_t& x) const;
    void            shift(colvec_t& x) const;

    mat_t           get_scale() const;
    colvec_t        get_shift() const;

    int m_nrows;
    real_t m_shift;
};

