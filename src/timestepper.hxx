#include "types.hxx"

class Timestepper {
public:
    Timestepper(colvec_t a_ic, real_t a_dt);
    void take_timestep(const colvec_t& a_rhs);
    const colvec_t& get_state() const;
    real_t get_dt() const;

private:
    colvec_t m_curr_state;
    real_t m_dt;
};
