#pragma once
#include <types.hxx>

template <typename Derived>
class AffineTransform {
public:
    AffineTransform() {};

    colvec_t                            apply(const colvec_t& x) const;
    void                                apply(colvec_t& x) const;

    colvec_t                            scale(const colvec_t& x) const;
    colvec_t                            shift(const colvec_t& x) const;

    void                                scale(colvec_t& x) const;
    void                                shift(colvec_t& x) const;

    template <typename EigenDerived>
    Eigen::EigenBase<EigenDerived>      get_scale() const;
    colvec_t                            get_shift() const;

private:
    Derived& derived() {return *static_cast<Derived*>(this);};
};

#include "affine_transform.txx"
