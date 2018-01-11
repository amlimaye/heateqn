TEST_DIR=tests
EXAMPLE_DIR=examples
SRC_DIR=src
INCL_DIR=include
BUILD_DIR=build

EIGEN_INCL_DIR=/usr/local/Cellar/eigen/3.3.4/include/eigen3

GTEST_INCL_DIR=/usr/local/include/
GTEST_LINK_DIR=/usr/local/lib/
GTEST_LDFLAGS=-L$(GTEST_LINK_DIR) -lgtest

PYBIND_INCL_DIR=$(shell python -m pybind11 --includes)
PYBIND_FLAGS=-shared -fPIC -undefined dynamic_lookup $(PYBIND_INCL_DIR)

CXX=g++
INCLUDES=-I$(EIGEN_INCL_DIR) -I$(INCL_DIR) -I$(GTEST_INCL_DIR)
CXXFLAGS=-std=c++11 -O3 -Wall -Wextra -Werror -g

COPY=cp -r

#class object files
constant_shift: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx

laplace: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx

scale: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx

forward_euler: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx

backward_euler: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx

#unit tests for classes
laplace_test: Makefile laplace
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GTEST_LDFLAGS) -o $(BUILD_DIR)/$@ \
		$(TEST_DIR)/$@.cxx \
		$(addsuffix .o, $(addprefix $(BUILD_DIR)/, $(filter-out Makefile, $^)))

constant_shift_test: Makefile constant_shift 
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GTEST_LDFLAGS) -o $(BUILD_DIR)/$@ \
		$(TEST_DIR)/$@.cxx \
		$(addsuffix .o, $(addprefix $(BUILD_DIR)/, $(filter-out Makefile, $^)))

scale_test: Makefile scale
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

#instructive examples
pyexample_convergence_rate: Makefile forward_euler backward_euler scale
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(PYBIND_FLAGS) -o $(BUILD_DIR)/$@.so \
		$(EXAMPLE_DIR)/$(subst pyexample_,,$@)/$@.cxx \
		$(addsuffix .o, $(addprefix $(BUILD_DIR)/, $(filter-out Makefile, $^)))
	$(COPY) $(BUILD_DIR)/$@.so $(EXAMPLE_DIR)/$(subst pyexample_,,$@)/$@.so
	$(COPY) $(BUILD_DIR)/$@.so.dSYM $(EXAMPLE_DIR)/$(subst pyexample_,,$@)/$@.so.dSYM
	python $(EXAMPLE_DIR)/$(subst pyexample_,,$@)/$(subst pyexample_,,$@).py

alltests: Makefile laplace_test forward_euler_test backward_euler_test \
			constant_shift_test scale_test


clean:
	rm -rf $(BUILD_DIR)/*.dSYM $(BUILD_DIR)/*_test $(BUILD_DIR)/*.so \
		$(BUILD_DIR)/*.o
	rm -rf $(EXAMPLE_DIR)/*/*.so $(EXAMPLE_DIR)/*/*.dSYM
