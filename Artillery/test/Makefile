GTEST_DIR = $(HOME)/gtest-1.7.0

# Where to find user code.
USER_DIR = ../src
TEST_DIR = .


# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include

# Flags passed to the C++ compiler.
CXXFLAGS += -std=c++0x -g -Wall -Wextra -pthread -Wno-missing-field-initializers

#SFML flags
LDFLAGS = -lsfml-graphics -lsfml-window -lsfml-system

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS = artillery_test

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

# House-keeping build targets.

all : $(TESTS)

clean :
	rm -f $(TESTS) gtest.a gtest_main.a *.o 

get-points :
	./artillery_test

run: artillery_test
	./artillery_test

#########################################################################
#
#                           ALA KOSKE
#########################################################################
# Builds gtest.a and gtest_main.a.

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^



#########################################################################
#
#                           ALA KOSKE
#########################################################################



# Builds a sample test.  A test should link with either gtest.a or
# gtest_main.a, depending on whether it defines its own main()
# function.


calc.o : $(USER_DIR)/calc.cpp $(USER_DIR)/calc.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/calc.cpp 

weapon.o : $(USER_DIR)/weapon.cpp $(USER_DIR)/weapon.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/weapon.cpp

player.o : $(USER_DIR)/player.cpp $(USER_DIR)/player.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/player.cpp

map.o : $(USER_DIR)/map.cpp $(USER_DIR)/map.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/map.cpp

game.o : $(USER_DIR)/game.cpp $(USER_DIR)/game.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/game.cpp

test_source.o : $(TEST_DIR)/test_source.cpp \
                     $(USER_DIR)/calc.hpp $(USER_DIR)/weapon.hpp $(USER_DIR)/player.hpp $(USER_DIR)/map.hpp $(USER_DIR)/game.hpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(TEST_DIR)/test_source.cpp

artillery_test : calc.o weapon.o map.o player.o game.o test_source.o  gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ $(LDFLAGS) -o $@
