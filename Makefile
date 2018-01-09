CXX=g++
INCLUDES=-I/usr/local/Cellar/eigen/3.3.4/include/eigen3
CXXFLAGS=-std=c++11 -O3 -Wall -Wextra -Werror -g
SRC_DIR=src
BUILD_DIR=build

GTEST_LDFLAGS=-lgtest

constant_shift: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx

laplace: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx

laplace_test: Makefile laplace
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GTEST_LDFLAGS) -o $(BUILD_DIR)/$@ \
		$(SRC_DIR)/$@.cxx \
		$(addsuffix .o, $(addprefix $(BUILD_DIR)/, $(filter-out Makefile, $^)))

forward_euler: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx

forward_euler_test: Makefile forward_euler laplace constant_shift
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GTEST_LDFLAGS) -o $(BUILD_DIR)/$@ \
		$(SRC_DIR)/$@.cxx \
		$(addsuffix .o, $(addprefix $(BUILD_DIR)/, $(filter-out Makefile, $^)))

backward_euler: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx

backward_euler_test: Makefile backward_euler laplace
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GTEST_LDFLAGS) -o $(BUILD_DIR)/$@ \
		$(SRC_DIR)/$@.cxx \
		$(addsuffix .o, $(addprefix $(BUILD_DIR)/, $(filter-out Makefile, $^)))

build: Makefile timestepper laplace
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(BUILD_DIR)/main  \
	$(addprefix $(BUILD_DIR)/, $(addsuffix .o, $(filter-out $<, $^))) $(SRC_DIR)/driver.cxx

clean:
	rm -rf $(BUILD_DIR)/main.dSYM $(BUILD_DIR)/main $(BUILD_DIR)/*.o
