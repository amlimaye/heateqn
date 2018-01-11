#include <operators/scale.hxx>

Scale::Scale(const int a_nrows, const real_t a_alpha) : 
    m_nrows(a_nrows),
    m_alpha(a_alpha)
{};

colvec_t Scale::apply(const colvec_t& x) const {
    return this->scale(x) + this->shift(x);
};

void Scale::apply(colvec_t& x) const {
    this->scale(x);
    this->shift(x);
};

colvec_t Scale::scale(const colvec_t& x) const {
    return x * m_alpha;
};

colvec_t Scale::shift(const colvec_t&) const {
    auto out = colvec_t(m_nrows);
    for (int k = 0; k < m_nrows; k++) {
        out(k) = 0;
    }
    return out;
};

void Scale::scale(colvec_t& x) const {
    x *= m_alpha;
};

void Scale::shift(colvec_t&) const {};

mat_t Scale::get_scale() const {
    auto out = mat_t();
    out.setIdentity(m_nrows, m_nrows);
    out *= m_alpha;
    return out;
};

colvec_t Scale::get_shift() const {
    auto out = colvec_t();
    out.setZero(m_nrows);
    return out;
};
