CXX=g++
INCLUDES=-I/usr/local/Cellar/eigen/3.3.4/include/eigen3
CXXFLAGS=-std=c++11 -O0 -Wall -Wextra -Werror -g
SRC_DIR=src
BUILD_DIR=build

GTEST_LDFLAGS=-lgtest


laplace: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx

laplacetest: Makefile laplace
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GTEST_LDFLAGS) -o $(BUILD_DIR)/$@ \
		$(SRC_DIR)/$@.cxx \
		$(addsuffix .o, $(addprefix $(BUILD_DIR)/, $(filter-out Makefile, $^)))

forward_euler: Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $(BUILD_DIR)/$@.o $(SRC_DIR)/$@.cxx

build: Makefile timestepper laplace
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(BUILD_DIR)/main  \
	$(addprefix $(BUILD_DIR)/, $(addsuffix .o, $(filter-out $<, $^))) $(SRC_DIR)/driver.cxx

clean:
	rm -rf $(BUILD_DIR)/main.dSYM $(BUILD_DIR)/main $(BUILD_DIR)/*.o
