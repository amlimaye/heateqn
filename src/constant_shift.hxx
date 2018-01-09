#include "types.hxx"
#include "affine_transform.hxx"

class ConstantShift : public AffineTransform {
public:
    ConstantShift(const real_t a_shift, const int a_nrows);

    colvec_t        apply(const colvec_t& x) const;
    void            apply(colvec_t& x) const;

    colvec_t        scale(const colvec_t& x) const;
    colvec_t        shift(const colvec_t& x) const;

    void            scale(colvec_t& x) const;
    void            shift(colvec_t& x) const;

    mat_t           get_scale() const;
    colvec_t        get_shift() const;

    real_t m_shift;
    int m_nrows;
};

