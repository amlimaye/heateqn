#include <gtest/gtest.h>
#include <domain.hxx>
#include <random>

TEST(Domain2DTest, InitializationAndGetters) {
    int nsamples = 100;
    std::mt19937 gen(0);
    auto dist = std::uniform_real_distribution<> (-10.0, 10.0);
    std::vector<int> npoints_list = {1, 10, 20, 50, 100, 300};

    for (const auto& npoints : npoints_list) {
        for (int k = 0; k < nsamples; k++) {
            auto value = dist(gen);
            auto d = Domain2D(npoints, value);

            auto u = d.get_interior_data();
            auto bu = d.get_boundary_data();
            auto npoints_per_dim = d.get_npoints_per_dim();

            for (int k = 0; k < npoints*npoints; k++) {
                EXPECT_DOUBLE_EQ(u[k], value);
            }
            
            for (int k = 0; k < 4*npoints; k++) {
                EXPECT_DOUBLE_EQ(bu[k], value);
            }

            EXPECT_EQ(npoints_per_dim, npoints);
        }
    }
}

typedef std::function<real_t(real_t, real_t)> Functor;
void reference_map_implementation(int npoints, const Functor& f,
                                  colvec_t& u, colvec_t& bu) {
    u.resize(npoints * npoints);
    bu.resize(npoints * 4);

    //convert x,y index pair into a flat, x-major indexing scheme
    auto get_gidx = [&] (const int x_idx, const int y_idx) -> int {
        return x_idx + (npoints * y_idx);
    };

    //fill in all points in the domain with the sampled plane
    for (int k = 0; k < npoints; k++) {
        double x = (double) (k+1) / (npoints+1);
        for (int l = 0; l < npoints; l++) {
            double y = (double) (l+1) / (npoints+1);
            u(get_gidx(k, l)) = f(x, y);
        }
    }

    //fill in all points in the boundary
    //left boundary
    for (int k = 0; k < npoints; k++) {
        double y = (double) (k+1) / (npoints+1);
        bu(k) = f(0.0, y);
    }

    //bottom boundary
    for (int k = 0; k < npoints; k++) {
        double x = (double) (k+1) / (npoints+1);
        bu(npoints + k) = f(x, 0.0);
    }

    //right boundary
    for (int k = 0; k < npoints; k++) {
        double y = (double) (k+1) / (npoints+1);
        bu(npoints*2 + k) = f(1.0, y);
    }

    //top boundary
    for (int k = 0; k < npoints; k++) {
        double x = (double) (k+1) / (npoints+1);
        bu(npoints*3 + k) = f(x, 1.0);
    }
}

TEST(Domain2DTest, MapParaboloid) {
    int nsamples = 100;
    std::srand(0);
    std::mt19937 gen(0);
    auto dist = std::uniform_real_distribution<> (-10.0, 10.0);
    std::vector<int> npoints_list = {10, 50, 100, 500};

    for (const auto& npoints : npoints_list) {
        for (int k = 0; k < nsamples; k++) {
            auto f_param = [] (real_t x, real_t y, 
                               real_t curvature, colvec_t& offset_vector) {
                auto x_off = offset_vector(0);
                auto y_off = offset_vector(1);

                auto z = (x - x_off) * (x - x_off) + (y - y_off) * (y - y_off);
                z *= (curvature / 2.0);
                return z;
            };

            colvec_t offset_vector(2);
            offset_vector.setRandom();
            real_t curvature = dist(gen);
            auto f = [&] (real_t x, real_t y) {
                return f_param(x, y, curvature, offset_vector);
            };
   
            //call the implementation in Domain2D
            auto d = Domain2D(npoints);
            d.map_over_domain(f);
            auto domain_u = d.get_interior_data();
            auto domain_bu = d.get_boundary_data();

            //call reference implementation
            colvec_t u, bu;
            reference_map_implementation(npoints, f, u, bu);

            //check that the answer is identical
            for (int k = 0; k < npoints * npoints; k++) {
                EXPECT_DOUBLE_EQ(domain_u(k), u(k));
            }
            for (int k = 0; k < 4 * npoints; k++) {
                EXPECT_DOUBLE_EQ(domain_bu(k), bu(k));
            }
        }
    }
}

TEST(Domain2DTest, MapPlane) {
    int nsamples = 100;
    std::vector<int> npoints_list = {10, 50, 100, 500};
    std::srand(0);

    //using these for both offset and normal vectors
    std::vector<real_t> x_off_list= {1.0, 0.5, 1.57, 5.0};
    std::vector<real_t> y_off_list= {-1.0, 0.7, 2.68, -5.0};

    for (const auto& npoints : npoints_list) {
        for (int k = 0; k < nsamples; k++) {
            auto f_param = [] (real_t x, real_t y, 
                               colvec_t& normal_vector, 
                               colvec_t& offset_vector) {
                    assert(normal_vector.size() == 3);
                    assert(offset_vector.size() == 3);

                    auto rz = (-1.0 * (x - offset_vector(0)) *
                               normal_vector(0) - (y - offset_vector(1)) *
                               normal_vector(1)) / normal_vector(2) +
                               offset_vector(2);
                    return rz;
            };

            colvec_t offset_vector(3);
            colvec_t normal_vector(3);
            offset_vector.setRandom();
            normal_vector.setRandom();

            auto f = [&] (real_t x, real_t y) {
                return f_param(x, y, normal_vector, offset_vector);
            };

            //call the implementation in Domain2D
            auto d = Domain2D(npoints);
            d.map_over_domain(f);
            auto domain_u = d.get_interior_data();
            auto domain_bu = d.get_boundary_data();

            //call reference implementation
            colvec_t u, bu;
            reference_map_implementation(npoints, f, u, bu);

            //check that the answer is identical
            for (int k = 0; k < npoints * npoints; k++) {
                EXPECT_DOUBLE_EQ(domain_u(k), u(k));
            }
            for (int k = 0; k < 4 * npoints; k++) {
                EXPECT_DOUBLE_EQ(domain_bu(k), bu(k));
            }
        }
    }
}

void reference_boundary_corrector_implementation(int npoints_per_dim,
                                                 const colvec_t& boundary_values,
                                                 colvec_t& boundary_corrector) {
    //convert x,y index pair into a flat, x-major indexing scheme
    auto get_gidx = [&] (const int x_idx, const int y_idx) -> int {
        return x_idx + (npoints_per_dim * y_idx);
    };

    boundary_corrector.resize(npoints_per_dim * npoints_per_dim);
    boundary_corrector.setZero(npoints_per_dim * npoints_per_dim);
    
    /* boundary_values is a flat array of four different boundary condition
     * vectors, the indexing scheme is illustrated by this diagram:
     *
     *         (4)
     *      ----------
     *     |          |
     *     |          |
     * (1) |          | (3)
     *     |          |
     *     |          |
     *      ----------
     *         (2)
     *
     * each edge of the boundary array must have npoints_per_dim elements, so 
     * in total we must have 4 * npoints_per_dim elements
     */
    assert(boundary_values.size() == 4 * npoints_per_dim);
    
    //left boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        boundary_corrector(get_gidx(0, k)) += 
            boundary_values(k);
    }

    //bottom boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        boundary_corrector(get_gidx(k, 0)) += 
            boundary_values(npoints_per_dim + k);
    }

    //right boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        boundary_corrector(get_gidx(npoints_per_dim - 1, k)) += 
            boundary_values(npoints_per_dim * 2 + k);
    }

    //top boundary
    for (int k = 0; k < npoints_per_dim; k++) {
        boundary_corrector(get_gidx(k, npoints_per_dim - 1)) += 
            boundary_values(npoints_per_dim * 3 + k);
    }
}

TEST(Domain2DTest, BoundaryCorrector) {
    int nsamples = 100;
    std::vector<int> npoints_list = {10, 50, 100, 500};

    std::mt19937 gen(0);
    std::uniform_real_distribution<> dist(-10.0, 10.0);

    for (const auto& npoints : npoints_list) {
        for (int k = 0; k < nsamples; k++) {
            auto value = dist(gen);

            auto d = Domain2D(npoints, value);
            auto bu = d.get_boundary_data();
            auto boundary_corrector = d.get_boundary_corrector();

            colvec_t ref_boundary_corrector;
            reference_boundary_corrector_implementation(npoints, bu, 
                                                        ref_boundary_corrector);
            for (int k = 0; k < npoints * npoints; k++) {
                EXPECT_DOUBLE_EQ(boundary_corrector(k), ref_boundary_corrector(k));
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
