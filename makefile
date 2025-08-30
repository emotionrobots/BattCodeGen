# Build with numeric Jacobians (default): `make`
# To use CasADi Jacobians from PyBaMM codegen (if present): `make USE_CASADI_JAC=1`


CC = gcc
CFLAGS = -O3 -std=c11 -ffast-math -Wall -Wextra -pedantic 
LDFLAGS = -lm -lc
SRC_SIM = simulate.c integrator.c linalg.c batt_model.c
OBJ_SIM = $(SRC_SIM:.c=.o)


.PHONY: all clean run test
all: simulate


simulate: $(OBJ_SIM)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)


run: simulate
	./simulate | tee output.csv


clean:
	rm -f *.o simulate output.csv
