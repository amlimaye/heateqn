#include <gtest/gtest.h>
#include <limits>
#include <vector>
#include <tuple>
#include <operators/scale.hxx>

void make_linear_ic(int npoints, const real_t& right_bc, colvec_t& u) {
    u.resize(npoints);
    auto drop_per_interval = (1.0 - right_bc) / (npoints + 1);
    for (int l = 0; l < npoints; l++) {
        u[l] = 1.0 - (l+1)*drop_per_interval;
    }
};


TEST(ScaleTest, ScalingCorrectness) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    std::vector<real_t> alpha_list = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            auto scale = Scale(npoints, alpha);
            
            colvec_t u_old, u;
            make_linear_ic(npoints, 0.5, u_old);
            u = scale.scale(static_cast<const colvec_t>(u_old));

            for (int k = 0; k < npoints; k++)
                EXPECT_DOUBLE_EQ(u(k), alpha * u_old(k));
        }
    }
}

TEST(ScaleTest, ShiftReturnCorrectness) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    std::vector<real_t> alpha_list = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            auto scale = Scale(npoints, alpha);
            
            colvec_t u_old, u;
            make_linear_ic(npoints, 0.5, u_old);
            u = scale.shift(static_cast<const colvec_t>(u_old));

            //returned shift should return the zero vector
            for (int k = 0; k < npoints; k++)
                EXPECT_DOUBLE_EQ(u(k), 0);
        }
    }
}

TEST(ScaleTest, ShiftInplaceCorrectness) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    std::vector<real_t> alpha_list = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            auto scale = Scale(npoints, alpha);
            
            colvec_t u_old, u;
            make_linear_ic(npoints, 0.5, u);
            u_old = u;
            scale.shift(u);

            //inplace shift should return the same vector again
            for (int k = 0; k < npoints; k++)
                EXPECT_DOUBLE_EQ(u(k), u_old(k));
        }
    }
}

TEST(ScaleTest, ApplyCorrectness) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    std::vector<real_t> alpha_list = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            auto scale = Scale(npoints, alpha);
            
            colvec_t u_old, u;
            make_linear_ic(npoints, 0.5, u_old);
            u = scale.apply(static_cast<const colvec_t>(u_old));

            for (int k = 0; k < npoints; k++)
                EXPECT_DOUBLE_EQ(u(k), u_old(k)*alpha);
        }
    }
}

TEST(ScaleTest, InplaceApplyEquivalentToReturnApply) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    std::vector<real_t> alpha_list = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            colvec_t u_apply, u_scale, u_shift;
            make_linear_ic(npoints, 0.50, u_apply);
            //make_linear_ic(npoints, 0.50, u_scale);
            //make_linear_ic(npoints, 0.50, u_shift);

            auto scale = Scale(npoints, alpha);

            auto returned_apply = scale.apply(static_cast<const colvec_t>(u_apply));
            //auto returned_scale = scale.scale(static_cast<const colvec_t>(u_scale));
            //auto returned_shift = scale.shift(static_cast<const colvec_t>(u_shift));
            scale.apply(u_apply);
            //scale.apply(u_scale);
            //scale.apply(u_shift);

            for (int k = 0; k < npoints; k++) {
                EXPECT_DOUBLE_EQ(returned_apply(k), u_apply(k));
                //EXPECT_DOUBLE_EQ(returned_scale(k), u_scale(k));
                //EXPECT_DOUBLE_EQ(returned_shift(k), u_shift(k));
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
