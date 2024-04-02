CC=g++

FLAGS=-Wall -Werror -Wextra -std=c++17 -g

CHECKFL=-lcheck -lpthread

GCFLAGS=-fprofile-arcs -ftest-coverage

SRC=$(wildcard *.cpp *.h)

TESTS=$(wildcard tests/*.cpp)

OBJ=$(SRC:.cpp=.o)

TEST_OBJ=$(TESTS:.cpp=.o)

ifeq ($(shell uname -s), Linux)
  CHECKFL += -lrt -lsubunit -lm
endif

all: clean s21_matrix_oop.a test

s21_matrix_oop.a: $(OBJ)
	ar cr $@ $?

%.o: %.cpp
	$(CC) $(FLAGS) -c  $? -o $@

clang:
	clang-format -style=Google -n $(SRC) $(TESTS)	
	clang-format -style=Google -i $(SRC) $(TESTS)

test: clean
	$(CC) $(FLAGS) $(SRC) $(TESTS) -lgtest -pthread -o test_cases $(CHECKFL) $(GCFLAGS)
	chmod +x ./test_cases
	./test_cases

clean:
	rm -rf *.a test_cases coverage* html* *.o *.gcda *.gcno

valgrind: test
	sudo valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -s ./test_cases

