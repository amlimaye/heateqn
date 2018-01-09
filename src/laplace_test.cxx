#include <gtest/gtest.h>
#include <limits>
#include <vector>
#include <tuple>
#include "laplace.hxx"

#define EXPECT_NEAR_DIGITS(x,y,d) EXPECT_NEAR(x, y, std::pow(10, -d))

void make_expo_ic(int npoints, const real_t& right_bc, colvec_t& u) {
    u.resize(npoints);
    double decay_rate = -1.0 * std::log(right_bc);

    for (int k = 0; k < npoints; k++) {
        double x = (double) (k+1) / (npoints+1);
        u(k) = std::exp(-1.0 * decay_rate * x);
    } 
};

TEST(LaplaceOperator1DTest, SpacingCalculation) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    for (const auto& npoints : npoints_list) {
        auto lap = LaplaceOperator1D(npoints, 0.0);
        double expected_dx = 1.0 / (npoints + 1);
        EXPECT_DOUBLE_EQ(lap.get_dx(), expected_dx);
    }
}

TEST(LaplaceOperator1DTest, ApplyInplaceEquivalentToApplyReturn) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    std::vector<real_t> right_bcs = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& npoints : npoints_list) {
        for (const auto& right_bc : right_bcs) {
            colvec_t u;
            make_expo_ic(npoints, right_bc, u);

            auto lap = LaplaceOperator1D(npoints, right_bc);

            auto returned_value = lap.apply(static_cast<const colvec_t>(u));
            lap.apply(u);

            for (int k = 0; k < npoints; k++) {
                EXPECT_DOUBLE_EQ(returned_value(k), u(k));
            }
        }
    }
}

TEST(LaplaceOperator1DTest, ShiftInplaceEquivalentToShiftReturn) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    std::vector<real_t> right_bcs = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& npoints : npoints_list) {
        for (const auto& right_bc : right_bcs) {
            colvec_t u;
            make_expo_ic(npoints, right_bc, u);

            auto lap = LaplaceOperator1D(npoints, right_bc);

            auto returned_value = lap.shift(static_cast<const colvec_t>(u));
            lap.shift(u);

            for (int k = 0; k < npoints; k++) {
                EXPECT_DOUBLE_EQ(returned_value(k), u(k));
            }
        }
    }
}

TEST(LaplaceOperator1DTest, ScaleInplaceEquivalentToScaleReturn) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    std::vector<real_t> right_bcs = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& npoints : npoints_list) {
        for (const auto& right_bc : right_bcs) {
            colvec_t u;
            make_expo_ic(npoints, right_bc, u);

            auto lap = LaplaceOperator1D(npoints, right_bc);

            auto returned_value = lap.scale(static_cast<const colvec_t>(u));
            lap.scale(u);

            for (int k = 0; k < npoints; k++) {
                EXPECT_DOUBLE_EQ(returned_value(k), u(k));
            }
        }
    }
}
TEST(LaplaceOperator1DTest, UnscaledBoundaryTerm) {
    int npoints = 100;
    std::vector<real_t> right_bcs = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (const auto& right_bc : right_bcs) {
        //construct the laplace operator and get its boundary term
        auto lap = LaplaceOperator1D(npoints, right_bc);
        auto bc = lap.get_boundary_term();

        for (int l = 0; l < npoints; l++) {
            // left boundary is scaled to 1, so should just be 1/(dx)**2
            if (l == 0) {
                EXPECT_DOUBLE_EQ(bc[l], 1.0);
            }
            // right boundary should be right_bc/(dx)**2
            else if (l == npoints-1) {
                EXPECT_DOUBLE_EQ(bc[l], right_bc);
            }
            // if not on the boundary, this should be exactly zero
            else {
                EXPECT_DOUBLE_EQ(bc[l], 0.0);
            }
        }
    }
}

TEST(LaplaceOperator1DTest, UnscaledLaplacian) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    real_t right_bc = 0.0;
    for (const auto& npoints: npoints_list) {
        //construct the laplace operator and get its laplacian matrix
        auto lap = LaplaceOperator1D(npoints, right_bc);
        auto lap_mat = lap.get_laplacian();

        for (int k = 0; k < npoints; k++) {
            //check diagonal term
            EXPECT_DOUBLE_EQ(lap_mat(k, k), -2.0);

            //check off-diagonal terms
            if (k != 0) {
                EXPECT_DOUBLE_EQ(lap_mat(k, k-1), 1.0);
                EXPECT_DOUBLE_EQ(lap_mat(k-1, k), 1.0);
            }

            //all other elements should be zero
            for (int l = 0; l < npoints; l++) {
                if (abs(l - k) > 1) {
                    EXPECT_DOUBLE_EQ(lap_mat(k, l), 0.0);
                }
            }
        }
    }
}

TEST(LaplaceOperator1DTest, InternalScaleFactor) {
    std::vector<int> npoints_list = {5, 10, 50, 100};
    std::vector<real_t> alpha_list = {1.0, 5.0, 10.0, 100.0};
    real_t right_bc = 0.5;

    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            auto lap = LaplaceOperator1D(npoints, right_bc, alpha);

            auto dx = lap.get_dx();
            auto scale_factor = alpha * std::pow(dx, -2);

            //check the laplacian
            auto lap_matrix_unscaled = lap.get_laplacian();
            auto lap_matrix_scaled = lap.get_scale();
            for (int k = 0; k < npoints; k++) {
                for (int l = 0; l < npoints; l++) {
                    EXPECT_DOUBLE_EQ(lap_matrix_scaled(k, l), 
                                     scale_factor * lap_matrix_unscaled(k, l));
                }
            }

            //check the boundary term
            auto boundary_unscaled = lap.get_boundary_term();
            auto boundary_scaled = lap.get_shift();
            for (int k = 0; k < npoints; k++) {
                EXPECT_DOUBLE_EQ(boundary_scaled(k), scale_factor * boundary_unscaled(k));
            }
        }
    }
}

TEST(LaplaceOperator1DTest, LaplacianOfAStraightLine) {
    /*
     * double precision arithmetic gives us 15 digits of precision, but since
     * we multiply by (1/npoints)**2 each time and then add up the results, we
     * lose two digits of additional precision for each extra digit in npoints
     * this accounts for the amount of precision loss here.
     */
    std::vector<std::tuple<int, int>> npoints_and_precision_list = 
        {{2, 15}, {20, 13}, {200, 11}, {2000, 9}};

    //makes a vector whose elements decay linearly on the interval (1.0, 0.0)
    auto mkstraightline = [](int npoints, colvec_t& u) {
        u.resize(npoints);
        auto drop_per_interval = 1.0 / (npoints + 1);
        for (int l = 0; l < npoints; l++) {
            u[l] = 1.0 - (l+1)*drop_per_interval;
        }
        return u;
    };

    for (const auto& elem : npoints_and_precision_list) {
        auto npoints = std::get<0>(elem);
        auto digits = std::get<1>(elem);
        auto lap = LaplaceOperator1D(npoints, 0.0);

        colvec_t u;
        mkstraightline(npoints, u);
        lap.apply(u);

        for (int k = 0; k < npoints; k++) {
            EXPECT_NEAR_DIGITS(u[k], 0.0, digits);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
