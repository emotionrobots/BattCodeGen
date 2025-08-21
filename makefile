BIN_DIR = bin
INC_DIR = include 
SRC_DIR = src
LIB_DIR = lib

CC = gcc
AR = ar
CFLAGS = -Wall -Wno-unused-parameter -Wno-sign-compare -Wextra -O3 -std=c11 -DNDEBUG -fPIC -I. -I$(INC_DIR)
LDFLAGS = 
LIBS = -lm -lc 


# Source files for binaries
APP_SRCS = $(SRC_DIR)/main.c 


# Source files for library
LIB_SRCS = $(SRC_DIR)/batt_model.c 


# Binary targets
BIN = $(BIN_DIR)/run_model


# Library file
LIB = $(LIB_DIR)/libbatt_model.a


# Object target
OBJ = $(LIB_SRCS:.c=.o)


# All targets
.PHONY: all clean
all: $(BIN) $(LIB)


# Build library
$(LIB): $(OBJ)
	$(AR) rcs $@ $^


# Build binary
$(BIN): $(APP_SRCS) $(LIB)
	$(CC) $(CFLAGS) -o $(BIN) $(APP_SRCS) $(LIB) $(LIBS) 


# Compile source files to object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@


# Clean target
clean:
	rm -f $(BIN) $(LIB) $(SRC_DIR)/*.o
	rm -f $(BIN) $(LIB)
	rm -f dfn_gittr.csv 

