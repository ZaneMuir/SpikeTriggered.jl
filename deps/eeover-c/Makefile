CC = gcc
CFLAGS  = -Wall -g -fPIC
LIBCFLAGS = -Wall -shared
LDFLAGS = -lgsl -lgslcblas

OBJ = solvers.o zsolve_quartic.o Roots3And4.o
ifeq ($(OS), Windows_NT)
	LIBNAME = libeeover.dll
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S), Linux)
		LIBNAME = libeeover.so
	endif
	ifeq ($(UNAME_S), Darwin)
		LIBNAME = libeeover.dylib
	endif
endif

library: $(OBJ)
	$(CC) $(OBJ) $(LIBCFLAGS) -o $(LIBNAME) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

julia-lib: library clean

.PHONY: clean
clean:
	rm -rf $(BIN) $(OBJ) *.gcda *gcno
clean-all: clean
	rm -rf $(LIBNAME)
