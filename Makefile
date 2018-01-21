TEST_DIR=tests
EXAMPLE_DIR=examples
SRC_DIR=src
INCL_DIR=include
BUILD_DIR=build
DEP_DIR=$(BUILD_DIR)/dep

VPATH = $(BUILD_DIR):$(SRC_DIR):$(INCL_DIR):$(TEST_DIR)

EIGEN_INCL_DIR=/usr/local/Cellar/eigen/3.3.4/include/eigen3

GTEST_INCL_DIR=/usr/local/include/
GTEST_LINK_DIR=/usr/local/lib/
GTEST_LDFLAGS=-L$(GTEST_LINK_DIR) -lgtest

PYBIND_INCL_DIR=$(shell python -m pybind11 --includes)
PYBIND_FLAGS=-shared -fPIC -undefined dynamic_lookup $(PYBIND_INCL_DIR)

CXX=g++
INCLUDES=-I$(EIGEN_INCL_DIR) -I$(INCL_DIR) -I$(GTEST_INCL_DIR)
CXXFLAGS=-std=c++14 -O3 -Wall -Wextra -Werror -g
DEPFLAGS=-MM

COPY=cp -r

#class object files
%.o: Makefile $(SRC_DIR)/%.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@ $(filter-out $<, $^)

#domain.hxx is templated and hence header-only. haven't figured out how to do
#the dependencies properly for this one
domain_utest: Makefile domain_test.cxx 
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GTEST_LDFLAGS) -o $(BUILD_DIR)/$@ \
		$(filter %.cxx, $^)

#unit tests for classes
%_utest: Makefile %_test.cxx %.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GTEST_LDFLAGS) -o $(BUILD_DIR)/$@ \
		$(filter %.cxx, $^) \
		$(addprefix $(BUILD_DIR)/, $(notdir $(filter %.o, $^)))

#some classes have unit tests that depend on other classes
TEST_OBJS=laplace.o constant_shift.o
%_test: Makefile %.o $(TEST_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GTEST_LDFLAGS) -o $(BUILD_DIR)/$@ \
		$(TEST_DIR)/$@.cxx \
		$(addprefix $(BUILD_DIR)/, $(filter-out Makefile, $^))

#instructive examples
OBJS=forward_euler.o backward_euler.o scale.o
pyexample_convergence_rate: Makefile $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(PYBIND_FLAGS) -o $(BUILD_DIR)/$@.so \
		$(EXAMPLE_DIR)/$(subst pyexample_,,$@)/$@.cxx \
		$(addprefix $(BUILD_DIR)/, $(filter-out Makefile, $^))
	$(COPY) $(BUILD_DIR)/$@.so $(EXAMPLE_DIR)/$(subst pyexample_,,$@)/$@.so
	$(COPY) $(BUILD_DIR)/$@.so.dSYM $(EXAMPLE_DIR)/$(subst pyexample_,,$@)/$@.so.dSYM
	python $(EXAMPLE_DIR)/$(subst pyexample_,,$@)/$(subst pyexample_,,$@).py

alltests: Makefile laplace_utest forward_euler_test backward_euler_test \
			constant_shift_utest scale_utest

allpyexamples: pyexample_convergence_rate

clean:
	rm -rf $(BUILD_DIR)/*.dSYM $(BUILD_DIR)/*_test $(BUILD_DIR)/*.so \
		$(BUILD_DIR)/*.o
	rm -rf $(EXAMPLE_DIR)/*/*.so $(EXAMPLE_DIR)/*/*.dSYM

.PHONY: pyexample_convergence_rate alltests allpyexamples clean
