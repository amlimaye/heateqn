#include <gtest/gtest.h>
#include <limits>
#include <vector>
#include <tuple>
#include "laplace.hxx"

TEST(LaplaceOperator1DTest, SpacingCalculation) {
    std::vector<int> npoints_list = {1, 2, 5, 10, 100, 1000};
    for (auto&& npoints : npoints_list) {
        auto lap = LaplaceOperator1D(npoints, 0.0);
        double expected_dx = 1.0 / (npoints + 1);
        EXPECT_DOUBLE_EQ(lap.get_dx(), expected_dx);
    }
}

TEST(LaplaceOperator1DTest, BoundaryTerm) {
    std::vector<int> npoints_list = {5, 10, 50, 100};
    std::vector<real_t> right_bcs = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (auto&& npoints : npoints_list) {
        for (auto&& right_bc : right_bcs) {
            //construct the laplace operator and get its boundary term
            auto lap = LaplaceOperator1D(npoints, right_bc);
            auto bc = lap.get_boundary_term();
            auto scale = std::pow(lap.get_dx(), -2);

            for (int l = 0; l < npoints; l++) {
                // left boundary is scaled to 1, so should just be 1/(dx)**2
                if (l == 0) {
                    EXPECT_DOUBLE_EQ(bc[l], scale);
                }
                // right boundary should be right_bc/(dx)**2
                else if (l == npoints-1) {
                    EXPECT_DOUBLE_EQ(bc[l], right_bc*scale);
                }
                // if not on the boundary, this should be exactly zero
                else {
                    EXPECT_DOUBLE_EQ(bc[l], 0.0);
                }
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
        {{2, 14}, {20, 12}, {200, 10}, {2000, 8}};

    //makes a vector whose elements decay linearly on the interval (1.0, 0.0)
    auto mkstraightline = [](int npoints, colvec_t& u) {
        u.resize(npoints);
        auto drop_per_interval = 1.0 / (npoints + 1);
        for (int l = 0; l < npoints; l++) {
            u[l] = 1.0 - (l+1)*drop_per_interval;
        }
        return u;
    };

    for (auto&& elem : npoints_and_precision_list) {
        auto npoints = std::get<0>(elem);
        auto digits = std::get<1>(elem);
        auto lap = LaplaceOperator1D(npoints, 0.0);

        colvec_t u;
        mkstraightline(npoints, u);
        auto Lu = lap.apply(u);

        for (int k = 0; k < npoints; k++) {
            EXPECT_NEAR(Lu[k], 0.0, std::pow(10, -digits));
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
