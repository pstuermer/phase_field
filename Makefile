CC = gcc

FFTW_INC = $(HOME)/fftw3/include
FFTW_LIB = $(HOME)/fftw3/lib

CFLAGS = -std=c11 -O3 -march=native -mtune=native -Wall -Wextra -Iinclude -funroll-loops -flto -fPIC #-fsanitize=address -g
LDFLAGS = -L./build -lfftw3 -L$(LAPACKE_LIB) -lm -flto #-fsanitize=address -g

#COPY_MAIN_SRC = $(wildcard *.c)

# does wildcard search for anything with .c ending in given directories
# and turns it into a list saved in SRC
SRC = $(wildcard src/*.c src/**/*.c)

# takes SRC list and substitues all .c files to .o to create list of
# objects to be created
OBJ = $(patsubst %.c, %.o,$(SRC))

//TEST_SRC = $(wildcard tests/apply_block_conj_grad_test.c)
TESTS = $(patsubst %.c,%.ex ,$(TEST_SRC))
TEST_OBJ = $(patsubst %.c, %.o, $(TEST_SRC))

TRG = build/libphase.a
SO_TRG = $(patsubst %.a, $.so, $(TRG))

all: $(TRG) $(SO_TRG) $(TESTS)

$(TRG): build $(OBJ)
	ar rcs $@ $(OBJ)
	ranlib $@

$(SO_TRG): $(TRG) $(OBJ)
	$(CC) -shared -o $@ $(OBJ)

build:
	@mkdir -p build
	@mkdir -p bin

# Build test executables
tests/%.o: tests/%.c
	$(CC) $(CFLAGS) -c $< -o $@

tests/%.ex: tests/%.o $(TRG)
	$(CC) $< -o $@ $(LDFLAGS) -fopenmp
	rm -f $<

# Run Valgrind on the specified test executable
valgrind: $(TESTS)
	@echo "Running Valgrind on all tests ..."
	@for test in $(TESTS); do \
		echo "Running Valgrind on $$test..."; \
		valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./$$test; \
	done

# The cleaner
clean:
	rm -rf build $(OBJ) $(TEST_OBJ) $(TESTS)

