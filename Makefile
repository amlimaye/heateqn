TEST_DIR=tests
SRC_DIR=src
INCL_DIR=include
BUILD_DIR=build

EIGEN_INCL_DIR=/usr/local/Cellar/eigen/3.3.4/include/eigen3

GTEST_INCL_DIR=/usr/local/include/
GTEST_LINK_DIR=/usr/local/lib/
GTEST_LDFLAGS=-L$(GTEST_LINK_DIR) -lgtest

CXX=g++
INCLUDES=-I$(EIGEN_INCL_DIR) -I$(INCL_DIR) -I$(GTEST_INCL_DIR)
CXXFLAGS=-std=c++11 -O3 -Wall -Wextra -Werror -g

constant_shift: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx

laplace: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx

forward_euler: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx

backward_euler: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx


laplace_test: Makefile laplace
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GTEST_LDFLAGS) -o $(BUILD_DIR)/$@ \
		$(TEST_DIR)/$@.cxx \
		$(addsuffix .o, $(addprefix $(BUILD_DIR)/, $(filter-out Makefile, $^)))

forward_euler_test: Makefile forward_euler laplace constant_shift
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GTEST_LDFLAGS) -o $(BUILD_DIR)/$@ \
		$(TEST_DIR)/$@.cxx \
		$(addsuffix .o, $(addprefix $(BUILD_DIR)/, $(filter-out Makefile, $^)))

backward_euler_test: Makefile backward_euler laplace constant_shift
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GTEST_LDFLAGS) -o $(BUILD_DIR)/$@ \
		$(TEST_DIR)/$@.cxx \
		$(addsuffix .o, $(addprefix $(BUILD_DIR)/, $(filter-out Makefile, $^)))

alltests: Makefile laplace_test forward_euler_test backward_euler_test

clean:
	rm -rf $(BUILD_DIR)/main.dSYM $(BUILD_DIR)/main $(BUILD_DIR)/*.o
