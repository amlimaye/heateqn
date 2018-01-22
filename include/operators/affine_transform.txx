template <typename Derived>
colvec_t AffineTransform<Derived>::apply(const colvec_t& x) const {
    derived().apply(x);
}

template <typename Derived>
void AffineTransform<Derived>::apply(colvec_t& x) const {
    derived().apply(x);
}

template <typename Derived>
colvec_t AffineTransform<Derived>::scale(const colvec_t& x) const {
    derived().scale(x);
}

template <typename Derived>
colvec_t AffineTransform<Derived>::shift(const colvec_t& x) const {
    derived().shift(x);
}

template <typename Derived>
void AffineTransform<Derived>::scale(colvec_t& x) const {
    derived().scale(x);
}

template <typename Derived>
void AffineTransform<Derived>::shift(colvec_t& x) const {
    derived().shift(x);
}

template <typename Derived>
template <typename EigenDerived>
Eigen::EigenBase<EigenDerived> AffineTransform<Derived>::get_scale() const {
    derived().get_scale();
};

template <typename Derived>
colvec_t AffineTransform<Derived>::get_shift() const {
    derived().get_scale();
};
