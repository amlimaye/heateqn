#include <gtest/gtest.h>
#include <limits>
#include <vector>
#include <tuple>
#include <operators/constant_shift.hxx>

void make_linear_ic(int npoints, const real_t& right_bc, colvec_t& u) {
    u.resize(npoints);
    auto drop_per_interval = (1.0 - right_bc) / (npoints + 1);
    for (int l = 0; l < npoints; l++) {
        u[l] = 1.0 - (l+1)*drop_per_interval;
    }
};


TEST(ConstantShiftTest, ScalingInplaceCorrectness) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    std::vector<real_t> alpha_list = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            auto op = ConstantShift(npoints, alpha);
            
            colvec_t u;
            make_linear_ic(npoints, 0.5, u);
            op.scale(u);

            //vector should get scaled to zero
            for (int k = 0; k < npoints; k++)
                EXPECT_DOUBLE_EQ(u(k), 0.0);
        }
    }
}

TEST(ConstantShiftTest, ScalingReturnCorrectness) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    std::vector<real_t> alpha_list = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            auto op = ConstantShift(npoints, alpha);
            
            colvec_t u, u_old;
            make_linear_ic(npoints, 0.5, u_old);
            u = op.scale(static_cast<const colvec_t>(u_old));

            //should return the zero vector
            for (int k = 0; k < npoints; k++)
                EXPECT_DOUBLE_EQ(u(k), 0.0);
        }
    }
}

TEST(ConstantShiftTest, ShiftReturnCorrectness) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    std::vector<real_t> alpha_list = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            auto op = ConstantShift(npoints, alpha);
            
            colvec_t u_old, u;
            make_linear_ic(npoints, 0.5, u_old);
            u = op.shift(static_cast<const colvec_t>(u_old));

            //returned shift should be the requested shift
            for (int k = 0; k < npoints; k++)
                EXPECT_DOUBLE_EQ(u(k), alpha);
        }
    }
}

TEST(ConstantShiftTest, ShiftInplaceCorrectness) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    std::vector<real_t> alpha_list = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            auto scale = ConstantShift(npoints, alpha);
            
            colvec_t u;
            make_linear_ic(npoints, 0.5, u);
            scale.shift(u);

            //inplace shift should be the requested shift
            for (int k = 0; k < npoints; k++)
                EXPECT_DOUBLE_EQ(u(k), alpha);
        }
    }
}

TEST(ConstantShiftTest, ApplyCorrectness) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    std::vector<real_t> alpha_list = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            auto op = ConstantShift(npoints, alpha);
            
            colvec_t u;
            make_linear_ic(npoints, 0.5, u);
            op.apply(u);

            //applied result should be the requested constant
            for (int k = 0; k < npoints; k++)
                EXPECT_DOUBLE_EQ(u(k), alpha);
        }
    }
}

TEST(ConstantShiftTest, InplaceApplyEquivalentToReturnApply) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    std::vector<real_t> alpha_list = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            colvec_t u_apply;
            make_linear_ic(npoints, 0.50, u_apply);

            auto op = ConstantShift(npoints, alpha);

            auto returned_apply = op.apply(static_cast<const colvec_t>(u_apply));
            op.apply(u_apply);

            for (int k = 0; k < npoints; k++) {
                EXPECT_DOUBLE_EQ(returned_apply(k), u_apply(k));
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
