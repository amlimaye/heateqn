#include <gtest/gtest.h>
#include <limits>
#include <vector>
#include <tuple>
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

void make_parabolic_ic_and_bc(int npoints_per_dim, 
                              const real_t curvature,
                              const colvec_t& offset_vector, 
                              colvec_t& u, colvec_t& bu) {
    //resize u and bu to the right sizes
    u.resize(npoints_per_dim * npoints_per_dim);
    bu.resize(npoints_per_dim * 4);

    //convert x,y index pair into a flat, x-major indexing scheme
    auto get_gidx = [&] (const int x_idx, const int y_idx) -> int {
        return x_idx + (npoints_per_dim * y_idx);
    };

    //fill in all points in the domain with the sampled plane
    for (int k = 0; k < npoints_per_dim; k++) {
        double x = (double) (k+1) / (npoints_per_dim+1);
        for (int l = 0; l < npoints_per_dim; l++) {
            double y = (double) (l+1) / (npoints_per_dim+1);
            u(get_gidx(k, l)) = sample_paraboloid(x, y, curvature,
                                                  offset_vector);
        }
    }

    //fill in all points in the boundary
    //left boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        double y = (double) (k+1) / (npoints_per_dim+1);
        bu(k) = sample_paraboloid(0.0, y, curvature, offset_vector);
    }

    //bottom boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        double x = (double) (k+1) / (npoints_per_dim+1);
        bu(npoints_per_dim + k) = 
            sample_paraboloid(x, 0.0, curvature, offset_vector);
    }

    //right boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        double y = (double) (k+1) / (npoints_per_dim+1);
        bu(npoints_per_dim*2 + k) = 
            sample_paraboloid(1.0, y, curvature, offset_vector);
    }

    //top boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        double x = (double) (k+1) / (npoints_per_dim+1);
        bu(npoints_per_dim*3 + k) = 
            sample_paraboloid(x, 1.0, curvature, offset_vector);
    }
}

void make_planar_ic_and_bc(int npoints_per_dim, 
                           const colvec_t& normal_vector,
                           const colvec_t& offset_vector,
                           colvec_t& u, colvec_t& bu) {
    //resize u and bu to the right sizes
    u.resize(npoints_per_dim * npoints_per_dim);
    bu.resize(npoints_per_dim * 4);

    //convert x,y index pair into a flat, x-major indexing scheme
    auto get_gidx = [&] (const int x_idx, const int y_idx) -> int {
        return x_idx + (npoints_per_dim * y_idx);
    };

    //fill in all points in the domain with the sampled plane
    for (int k = 0; k < npoints_per_dim; k++) {
        double x = (double) (k+1) / (npoints_per_dim+1);
        for (int l = 0; l < npoints_per_dim; l++) {
            double y = (double) (l+1) / (npoints_per_dim+1);
            u(get_gidx(k, l)) = sample_plane(x, y, normal_vector, offset_vector);
        }
    }

    //fill in all points in the boundary with the sampled plane
    
    //left boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        double y = (double) (k+1) / (npoints_per_dim+1);
        bu(k) = sample_plane(0.0, y, normal_vector, offset_vector);
    }

    //bottom boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        double x = (double) (k+1) / (npoints_per_dim+1);
        bu(npoints_per_dim + k) = 
            sample_plane(x, 0.0, normal_vector, offset_vector);
    }

    //right boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        double y = (double) (k+1) / (npoints_per_dim+1);
        bu(npoints_per_dim*2 + k) = 
            sample_plane(1.0, y, normal_vector, offset_vector);
    }

    //top boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        double x = (double) (k+1) / (npoints_per_dim+1);
        bu(npoints_per_dim*3 + k) = 
            sample_plane(x, 1.0, normal_vector, offset_vector);
    }
}

TEST(LaplaceOperator2DTest, SpacingCalculation) {
    std::vector<int> npoints_list = {5};

    colvec_t normal_vector;
    normal_vector.resize(3);
    normal_vector << 0.0, 1.0, 1.0;
    colvec_t offset_vector;
    offset_vector.resize(3);
    offset_vector << 0.0, 0.0, 1.0;

    for (const auto& npoints : npoints_list) {
        colvec_t u, bu;
        make_planar_ic_and_bc(npoints, normal_vector, offset_vector, u, bu);

        auto lap = LaplaceOperator2D(npoints, bu);
        double expected_dx = 1.0 / (npoints + 1);
        EXPECT_DOUBLE_EQ(lap.get_dx(), expected_dx);
    }
}

TEST(LaplaceOperator2DTest, UnscaledLaplacian) {
    std::vector<int> npoints_list = {5, 10, 50, 100};

    colvec_t normal_vector;
    normal_vector.resize(3);
    normal_vector << 0.0, 1.0, 1.0;
    colvec_t offset_vector;
    offset_vector.resize(3);
    offset_vector << 0.0, 0.0, 1.0;

    for (const auto& npoints: npoints_list) {
        colvec_t u, bu; 
        make_planar_ic_and_bc(npoints, normal_vector, offset_vector, u, bu);

        //construct the laplace operator and get its laplacian matrix
        auto lap = LaplaceOperator2D(npoints, bu);
        auto lap_mat = lap.get_laplacian();

        //reference implementation of the algorithm
        //convert x,y index pair into a flat, x-major indexing scheme
        auto get_gidx = [&] (const int x_idx, const int y_idx) -> int {
            return x_idx + (npoints * y_idx);
        };

        mat_t reference_matrix;
        reference_matrix.resize(npoints*npoints, npoints*npoints);
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

        for (int k = 0; k < npoints*npoints; k++) {
            for (int l = 0; l < npoints*npoints; l++) {
                EXPECT_DOUBLE_EQ(lap_mat(k, l), reference_matrix(k, l));
            }
        }
    }
}

TEST(LaplaceOperator2DTest, UnscaledBoundaryTerm) {
    std::vector<int> npoints_list = {5, 10, 50, 100};

    colvec_t normal_vector;
    normal_vector.resize(3);
    normal_vector << 0.0, 1.0, 1.0;
    colvec_t offset_vector;
    offset_vector.resize(3);
    offset_vector << 0.0, 0.0, 1.0;

    for (const auto& npoints: npoints_list) {
        colvec_t u, bu; 
        make_planar_ic_and_bc(npoints, normal_vector, offset_vector, u, bu);

        //construct the laplace operator and get its laplacian matrix
        auto lap = LaplaceOperator2D(npoints, bu);
        auto boundary_term = lap.get_boundary_term();

        //reference implementation
        //convert x,y index pair into a flat, x-major indexing scheme
        auto get_gidx = [&] (const int x_idx, const int y_idx) -> int {
            return x_idx + (npoints * y_idx);
        };

        colvec_t ref_boundary_term;
        ref_boundary_term.resize(npoints*npoints);
        ref_boundary_term.setZero(npoints*npoints);
        
        assert(bu.size() == 4 * npoints);
        
        //left boundary
        for (int k = 0; k < npoints; k++) {
            ref_boundary_term(get_gidx(0, k)) += 
                bu(k);
        }
    
        //bottom boundary
        for (int k = 0; k < npoints; k++) {
            ref_boundary_term(get_gidx(k, 0)) += 
                bu(npoints + k);
        }
    
        //right boundary
        for (int k = 0; k < npoints; k++) {
            ref_boundary_term(get_gidx(npoints - 1, k)) += 
                bu(npoints * 2 + k);
        }
    
        //top boundary
        for (int k = 0; k < npoints; k++) {
            ref_boundary_term(get_gidx(k, npoints - 1)) += 
                bu(npoints * 3 + k);
        } 
        
        //check accuracy against reference
        for (int k = 0; k < npoints*npoints; k++) {
            EXPECT_DOUBLE_EQ(boundary_term(k), ref_boundary_term(k));
        }
    }
}

TEST(LaplaceOperator2DTest, LaplacianOfAParaboloid) {
    int npoints = 20;
    int digits = 9;

    std::vector<real_t> curvature_list = {-2.0, -1.0, 0.5, 1.0, 2.0, 5.0, 10.0};
    std::vector<real_t> x_offsets = {-5.0, -2.0, 0.0, 2.0, 5.0};
    std::vector<real_t> y_offsets = {-4.0, -1.5, 0.0, 1.0, 3.6};

    for (const auto& curvature : curvature_list) {
        for (const auto& x_o : x_offsets) {
            for (const auto& y_o : y_offsets) {
                colvec_t offset_vector;
                offset_vector.resize(2);
                offset_vector << x_o, y_o;
                
                colvec_t u, bu; 
                make_parabolic_ic_and_bc(npoints, curvature, 
                                         offset_vector, u, bu);
                
                auto lap = LaplaceOperator2D(npoints, bu);
                lap.apply(u);
 
                for (int k = 0; k < npoints*npoints; k++) {
                    EXPECT_NEAR_DIGITS(u[k], 2.0 * curvature, digits);
                }
            }
        }
    }
}

TEST(LaplaceOperator2DTest, LaplacianOfAPlane) {
    /* TODO: why am I losing so many digits here? in the 1D case, we had 13 
     * digits of precision for npoints = 20, but here we are doing 20x more 
     * adds, which cuts our precision by two more digits, but i see numerically 
     * that it's actually being cut by 4 digits! investigate.
     */
    int npoints = 20;
    int digits = 9;

    std::vector<real_t> x_normals = {-1.0, 1.2, 4.3, 5.6};
    std::vector<real_t> y_normals = {1.1, 1.3, -4.8, 5.2};
    std::vector<real_t> z_normals = {0.8, -1.5, 4.4, 5.1} ;

    std::vector<real_t> x_offsets = {1.0, -2.0, 5.0};
    std::vector<real_t> y_offsets = {1.1, 2.2, -5.6};
    std::vector<real_t> z_offsets = {-0.9, 2.1, 5.8};

    for (const auto& x_n : x_normals) {
        for (const auto& y_n : y_normals) {
            for (const auto& z_n : z_normals) {
                for (const auto& x_o : x_offsets) {
                    for (const auto& y_o : y_offsets) {
                        for (const auto& z_o : z_offsets) {
                            colvec_t normal_vector, offset_vector;
                            normal_vector.resize(3);
                            offset_vector.resize(3);
                            
                            normal_vector << x_n, y_n, z_n;
                            offset_vector << x_o, y_o, z_o;

                            colvec_t u, bu; 
                            make_planar_ic_and_bc(npoints, normal_vector, 
                                    offset_vector, u, bu);
                    
                            auto lap = LaplaceOperator2D(npoints, bu);
                            lap.apply(u);
                    
                            for (int k = 0; k < npoints*npoints; k++) {
                                EXPECT_NEAR_DIGITS(u[k], 0.0, digits);
                            }
                        }
                    }
                }
            }
        }
    }
}

TEST(LaplaceOperator2DTest, ApplyInplaceEquivalentToApplyReturn) {
    std::vector<int> npoints_list = {3, 5, 10, 20, 50};
    colvec_t normal_vector;
    normal_vector.resize(3);
    normal_vector << 0.0, 1.0, 1.0;
    colvec_t offset_vector;
    offset_vector.resize(3);
    offset_vector << 0.0, 0.0, 1.0;
    
    for (const auto& npoints : npoints_list) {
        colvec_t u, bu; 
        make_planar_ic_and_bc(npoints, normal_vector, offset_vector, u, bu);

        auto lap = LaplaceOperator2D(npoints, bu);

        auto returned_value = lap.apply(static_cast<const colvec_t>(u));
        lap.apply(u);

        for (int k = 0; k < npoints*npoints; k++) {
            EXPECT_DOUBLE_EQ(returned_value(k), u(k));
        }
    }
}

TEST(LaplaceOperator2DTest, ShiftInplaceEquivalentToShiftReturn) {
    std::vector<int> npoints_list = {3, 5, 10, 20, 50};
    colvec_t normal_vector;
    normal_vector.resize(3);
    normal_vector << 0.0, 1.0, 1.0;
    colvec_t offset_vector;
    offset_vector.resize(3);
    offset_vector << 0.0, 0.0, 1.0;
    
    for (const auto& npoints : npoints_list) {
        colvec_t u, bu; 
        make_planar_ic_and_bc(npoints, normal_vector, offset_vector, u, bu);

        auto lap = LaplaceOperator2D(npoints, bu);

        auto returned_value = lap.shift(static_cast<const colvec_t>(u));
        lap.shift(u);

        for (int k = 0; k < npoints*npoints; k++) {
            EXPECT_DOUBLE_EQ(returned_value(k), u(k));
        }
    }
}

TEST(LaplaceOperator2DTest, ScaleInplaceEquivalentToScaleReturn) {
    std::vector<int> npoints_list = {3, 5, 10, 20, 50};
    colvec_t normal_vector;
    normal_vector.resize(3);
    normal_vector << 0.0, 1.0, 1.0;
    colvec_t offset_vector;
    offset_vector.resize(3);
    offset_vector << 0.0, 0.0, 1.0;
    
    for (const auto& npoints : npoints_list) {
        colvec_t u, bu; 
        make_planar_ic_and_bc(npoints, normal_vector, offset_vector, u, bu);

        auto lap = LaplaceOperator2D(npoints, bu);

        auto returned_value = lap.scale(static_cast<const colvec_t>(u));
        lap.scale(u);

        for (int k = 0; k < npoints*npoints; k++) {
            EXPECT_DOUBLE_EQ(returned_value(k), u(k));
        }
    }
}

TEST(LaplaceOperator2DTest, InternalScaleFactor) {
    std::vector<int> npoints_list = {3, 5, 10, 20, 50};
    std::vector<real_t> alpha_list = {1.0, 5.0, 10.0, 100.0};

    colvec_t normal_vector;
    normal_vector.resize(3);
    normal_vector << 0.0, 1.0, 1.0;
    colvec_t offset_vector;
    offset_vector.resize(3);
    offset_vector << 0.0, 0.0, 1.0;
 
    for (const auto& npoints : npoints_list) {
        for (const auto& alpha : alpha_list) {
            colvec_t u, bu; 
            make_planar_ic_and_bc(npoints, normal_vector, offset_vector, u, bu);
            auto lap = LaplaceOperator2D(npoints, bu, alpha);

            auto dx = lap.get_dx();
            auto scale_factor = alpha * std::pow(dx, -2);

            //check the laplacian
            auto lap_matrix_unscaled = lap.get_laplacian();
            auto lap_matrix_scaled = lap.get_scale();
            for (int k = 0; k < npoints*npoints; k++) {
                for (int l = 0; l < npoints*npoints; l++) {
                    EXPECT_DOUBLE_EQ(lap_matrix_scaled(k, l), 
                                     scale_factor * lap_matrix_unscaled(k, l));
                }
            }

            //check the boundary term
            auto boundary_unscaled = lap.get_boundary_term();
            auto boundary_scaled = lap.get_shift();
            for (int k = 0; k < npoints*npoints; k++) {
                EXPECT_DOUBLE_EQ(boundary_scaled(k), scale_factor * boundary_unscaled(k));
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
