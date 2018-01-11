#include <operators/constant_shift.hxx>

ConstantShift::ConstantShift(const int a_nrows, const real_t a_shift) : 
    m_nrows(a_nrows),
    m_shift(a_shift)
{};

colvec_t ConstantShift::apply(const colvec_t& x) const {
    return this->scale(x) + this->shift(x);
};

void ConstantShift::apply(colvec_t& x) const {
    this->scale(x);
    this->shift(x);
};

colvec_t ConstantShift::scale(const colvec_t& x) const {
    return x * 0.0;
};

colvec_t ConstantShift::shift(const colvec_t& x) const {
    auto out = colvec_t();
    out.setOnes(x.rows());
    out *= m_shift;
    return out;
};

void ConstantShift::scale(colvec_t& x) const {
    x *= 0.0;
};

void ConstantShift::shift(colvec_t& x) const {
    for (int k = 0; k < x.rows(); k++)
        x(k) = m_shift;
};

mat_t ConstantShift::get_scale() const {
    auto out = mat_t();
    out.setZero(m_nrows, m_nrows);
    return out;
};

colvec_t ConstantShift::get_shift() const {
    auto out = colvec_t();
    out.setOnes(m_nrows);
    out *= m_shift;
    return out;
};
