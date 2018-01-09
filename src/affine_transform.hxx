#pragma once
#include "types.hxx"

class AffineTransform {
public:
    AffineTransform() {};

    virtual colvec_t        apply(const colvec_t& x) const = 0;
    virtual void            apply(colvec_t& x) const = 0;

    virtual colvec_t        scale(const colvec_t& x) const = 0;
    virtual colvec_t        shift(const colvec_t& x) const = 0;

    virtual void            scale(colvec_t& x) const = 0;
    virtual void            shift(colvec_t& x) const = 0;
 
    virtual mat_t           get_scale() const = 0;
    virtual colvec_t        get_shift() const = 0;
};
