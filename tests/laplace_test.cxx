#include <gtest/gtest.h>
#include <limits>
#include <vector>
#include <tuple>
#include <random>
#include <operators/laplace.hxx>

#define EXPECT_NEAR_DIGITS(x,y,d) EXPECT_NEAR(x, y, std::pow(10, -d))

void make_expo_ic(int npoints, const real_t& right_bc, colvec_t& u) {
    u.resize(npoints);
    double decay_rate = -1.0 * std::log(right_bc);

    for (int k = 0; k < npoints; k++) {
        double x = (double) (k+1) / (npoints+1);
        u(k) = std::exp(-1.0 * decay_rate * x);
    } 
};

real_t sample_plane(const real_t x,
                    const real_t y,
                    const colvec_t& normal_vector,
                    const colvec_t& offset_vector) {
    assert(normal_vector.size() == 3);
    assert(offset_vector.size() == 3);

    auto rz = (-1.0 * (x - offset_vector(0)) * normal_vector(0)
               - (y - offset_vector(1)) * normal_vector(1)) / normal_vector(2)
              + offset_vector(2);
    return rz;
}

real_t sample_paraboloid(const real_t x,
                         const real_t y,
                         const real_t& curvature,
                         const colvec_t& offset_vector) {
    assert(offset_vector.size() == 2);
    auto x_off = offset_vector(0);
    auto y_off = offset_vector(1);

    auto z = (x - x_off) * (x - x_off) + (y - y_off) * (y - y_off);
    z *= (curvature / 2.0);
    return z;
}

class LaplaceOperator2DTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        std::srand(0);
        gen = std::mt19937(0);
        dist = std::uniform_real_distribution<>(-10.0, 10.0);

        normal_vector.resize(3);
        normal_vector << 0.0, 1.0, 1.0;

        offset_vector.resize(3);
        offset_vector << 0.0, 0.0, 1.0;
    }

    std::mt19937 gen;
    std::uniform_real_distribution<> dist;
    colvec_t normal_vector;
    colvec_t offset_vector;
};

TEST_F(LaplaceOperator2DTest, SpacingCalculation) {
    std::vector<int> npoints_list = {5, 10, 50, 100};

    for (const auto& npoints : npoints_list) {
        auto domain = Domain2D(npoints);
        auto lap = LaplaceOperator2D(domain);

        double expected_dx = 1.0 / (npoints + 1);
        EXPECT_DOUBLE_EQ(lap.get_dx(), expected_dx);
    }
}

void reference_laplacian_implementation(int npoints, mat_t& reference_matrix) {
    //convert x,y index pair into a flat, x-major indexing scheme
    auto get_gidx = [&] (const int x_idx, const int y_idx) -> int {
        return x_idx + (npoints * y_idx);
    };

    reference_matrix.resize(npoints*npoints, npoints*npoints);
    reference_matrix.setZero();
    for (int x_idx = 0; x_idx < npoints; x_idx++) {
        for (int y_idx = 0; y_idx < npoints; y_idx++) {
            //set the diagonal term
            auto diag_idx = get_gidx(x_idx, y_idx);
            reference_matrix(diag_idx, diag_idx) = -4;

            //set the term to the right
            if (x_idx < npoints-1) {
                auto right_idx = get_gidx(x_idx + 1, y_idx);
                reference_matrix(diag_idx, right_idx) = 1;
            }

            //set the term below
            if (y_idx < npoints-1) {
                auto down_idx = get_gidx(x_idx, y_idx + 1);
                reference_matrix(diag_idx, down_idx) = 1;
            }

            //set the term above
            if (y_idx > 0) {
                auto above_idx = get_gidx(x_idx, y_idx - 1);
                reference_matrix(diag_idx, above_idx) = 1;
            }

            //set the term to the left
            if (x_idx > 0) {
                auto left_idx = get_gidx(x_idx - 1, y_idx);
                reference_matrix(diag_idx, left_idx) = 1;
            }
        }
    }
}

TEST_F(LaplaceOperator2DTest, UnscaledLaplacian) {
    std::vector<int> npoints_list = {5, 10, 50};

    for (const auto& npoints: npoints_list) {
        auto domain = Domain2D(npoints);

        //construct the laplace operator and get its laplacian matrix
        auto lap = LaplaceOperator2D(domain);
        //explicitly cast to a dense matrix here because the ref impl. is dense
        auto lap_mat = mat_t(lap.get_laplacian());
        
        //call the reference implementation of the laplacian matrix
        mat_t ref_lap_mat;
        reference_laplacian_implementation(npoints, ref_lap_mat);

        //reference implementation of the algorithm
        for (int k = 0; k < npoints*npoints; k++) {
            for (int l = 0; l < npoints*npoints; l++) {
                EXPECT_DOUBLE_EQ(lap_mat(k, l), ref_lap_mat(k, l));
            }
        }
    }
}

TEST_F(LaplaceOperator2DTest, UnscaledBoundaryTerm) {
    std::vector<int> npoints_list = {5, 10, 50, 100};

    auto f = [&] (real_t x, real_t y) -> real_t {
        return sample_plane(x, y, normal_vector, offset_vector);
    };

    for (const auto& npoints: npoints_list) {
        //construct the domain
        auto domain = Domain2D(npoints);
        domain.map_over_domain(f);
        auto dbc = domain.get_boundary_corrector();

        //construct the laplace operator and get its laplacian matrix
        auto lap = LaplaceOperator2D(domain);
        auto boundary_term = lap.get_boundary_term();

        //check accuracy against reference
        for (int k = 0; k < npoints*npoints; k++) {
            EXPECT_DOUBLE_EQ(boundary_term(k), dbc(k));
        }
    }
}

TEST_F(LaplaceOperator2DTest, LaplacianOfAParaboloid) {
    int nsamples = 10;
    int npoints = 1000;
    int digits = 7;

    for (int k = 0; k < nsamples; k++) {
        //sample the paraboloid offset and curvature
        colvec_t offset_vector;
        offset_vector.resize(2);
        offset_vector.setRandom();
        real_t curvature = dist(gen);
      
        //prepare the domain 
        auto domain = Domain2D(npoints);
        auto f = [&] (real_t x, real_t y) {
            return sample_paraboloid(x, y, curvature, offset_vector);
        };
        domain.map_over_domain(f);

        //apply to the domain's interior data
        auto u = domain.get_interior_data();
        auto lap = LaplaceOperator2D(domain);
        lap.apply(u);
 
        for (int k = 0; k < npoints*npoints; k++) {
            //we expect 2*curvature because we add u_xx + u_yy
            EXPECT_NEAR_DIGITS(u[k], 2.0 * curvature, digits);
        }
    }
}

TEST_F(LaplaceOperator2DTest, LaplacianOfAPlane) {
    /* TODO: why am I losing so many digits here? in the 1D case, we had 13 
     * digits of precision for npoints = 20, but here we are doing 20x more 
     * adds, which cuts our precision by two more digits, but i see numerically 
     * that it's actually being cut by 4 digits! investigate.
     */
    int npoints = 1000;
    int digits = 7;
    int nsamples = 10;

    for (int k = 0; k < nsamples; k++) {
        //prepare the domain
        auto domain = Domain2D(npoints);

        //sample a plane over the domain
        colvec_t normal_vector(3);
        normal_vector.setRandom();
        colvec_t offset_vector(3);
        offset_vector.setRandom();
        auto f = [&] (real_t x, real_t y) {
            return sample_plane(x, y, normal_vector, offset_vector);
        };
        domain.map_over_domain(f);

        //pull out the interior data and apply the laplacian
        auto u = domain.get_interior_data();
        auto lap = LaplaceOperator2D(domain);
        lap.apply(u);
        
        for (int k = 0; k < npoints*npoints; k++) {
            EXPECT_NEAR_DIGITS(u[k], 0.0, digits);
        }
    }
}

TEST_F(LaplaceOperator2DTest, ApplyInplaceEquivalentToApplyReturn) {
    std::vector<int> npoints_list = {3, 5, 10, 20, 50};
   
    for (const auto& npoints : npoints_list) {
        //prepare the domain
        auto domain = Domain2D(npoints);
        auto f = [&] (real_t x, real_t y) {
            return sample_plane(x, y, normal_vector, offset_vector);
        };
        domain.map_over_domain(f);

        auto u = domain.get_interior_data();
        auto lap = LaplaceOperator2D(domain);
        auto returned_value = lap.apply(static_cast<const colvec_t>(u));
        lap.apply(u);

        for (int k = 0; k < npoints*npoints; k++) {
            EXPECT_DOUBLE_EQ(returned_value(k), u(k));
        }
    }
}

TEST_F(LaplaceOperator2DTest, ShiftInplaceEquivalentToShiftReturn) {
    std::vector<int> npoints_list = {3, 5, 10, 20, 50};
  
    for (const auto& npoints : npoints_list) {
        //prepare the domain
        auto domain = Domain2D(npoints);
        auto f = [&] (real_t x, real_t y) {
            return sample_plane(x, y, normal_vector, offset_vector);
        };
        domain.map_over_domain(f);

        auto u = domain.get_interior_data();
        auto lap = LaplaceOperator2D(domain);
        auto returned_value = lap.shift(static_cast<const colvec_t>(u));
        lap.shift(u);

        for (int k = 0; k < npoints*npoints; k++) {
            EXPECT_DOUBLE_EQ(returned_value(k), u(k));
        }
    }
}

TEST_F(LaplaceOperator2DTest, ScaleInplaceEquivalentToScaleReturn) {
    std::vector<int> npoints_list = {3, 5, 10, 20, 50};
   
    for (const auto& npoints : npoints_list) {
        //prepare the domain
        auto domain = Domain2D(npoints);
        auto f = [&] (real_t x, real_t y) {
            return sample_plane(x, y, normal_vector, offset_vector);
        };
        domain.map_over_domain(f);

        auto u = domain.get_interior_data();
        auto lap = LaplaceOperator2D(domain);
        auto returned_value = lap.scale(static_cast<const colvec_t>(u));
        lap.scale(u);

        for (int k = 0; k < npoints*npoints; k++) {
            EXPECT_DOUBLE_EQ(returned_value(k), u(k));
        }
    }
}

TEST_F(LaplaceOperator2DTest, InternalScaleFactor) {
    std::vector<int> npoints_list = {3, 5, 10, 20, 50};
    std::mt19937 gen(0);
    std::uniform_real_distribution<> dist (1.0, 100.0);
    int nsamples = 100;

    for (const auto& npoints : npoints_list) {
        for (int k = 0; k < nsamples; k++) {
            auto alpha = dist(gen);
            //prepare the domain
            auto domain = Domain2D(npoints);
            auto f = [&] (real_t x, real_t y) {
                return sample_plane(x, y, normal_vector, offset_vector);
            };
            domain.map_over_domain(f);

            auto u = domain.get_interior_data();
            auto lap = LaplaceOperator2D(domain, alpha);

            auto dx = lap.get_dx();
            auto scale_factor = alpha * std::pow(dx, -2);

            //check the laplacian
            auto lap_matrix_unscaled = lap.get_laplacian();
            auto lap_matrix_scaled = lap.get_scale();
    
            EXPECT_EQ(lap_matrix_unscaled.outerSize(), 
                      lap_matrix_scaled.outerSize());

            for (int k = 0; k < lap_matrix_scaled.outerSize(); k++) {
                {
                    sparse_mat_t::InnerIterator it_scaled(lap_matrix_scaled, k);
                    sparse_mat_t::InnerIterator it_unscaled(lap_matrix_unscaled, k);
                    for (; it_scaled; ++it_scaled, ++it_unscaled) {
                        EXPECT_DOUBLE_EQ(it_scaled.value(), 
                                         scale_factor*it_unscaled.value());
                    }
                }
            }

            /*
            for (int k = 0; k < npoints*npoints; k++) {
                for (int l = 0; l < npoints*npoints; l++) {
                    EXPECT_DOUBLE_EQ(lap_matrix_scaled(k, l), 
                                     scale_factor * lap_matrix_unscaled(k, l));
                }
            }*/

            //check the boundary term
            auto boundary_unscaled = lap.get_boundary_term();
            auto boundary_scaled = lap.get_shift();
            for (int k = 0; k < npoints*npoints; k++) {
                EXPECT_DOUBLE_EQ(boundary_scaled(k), 
                                 scale_factor * boundary_unscaled(k));
            }
        }
    }
}

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
    std::vector<int> npoints_list = {1, 2, 5, 10, 100};
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
    std::mt19937 gen(0);
    std::uniform_real_distribution<> dist (1.0, 100.0);
    int nsamples = 100;

    real_t right_bc = 0.5;

    for (const auto& npoints : npoints_list) {
        for (int k = 0; k < nsamples; k++) {
            auto alpha = dist(gen);
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
