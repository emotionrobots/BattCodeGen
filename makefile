BIN_DIR = bin
INC_DIR = include 
SRC_DIR = src

CC = gcc
AR = ar rcs
CFLAGS = -Wall -Wno-unused-parameter -Wno-sign-compare -Wextra -O3 -std=c11 -DNDEBUG -fPIC -I. -I$(INC_DIR)
LDFLAGS = 
LIBS = -lm -lc 
#LIBS = -lm -lc -lsundials_ida -lsundials_nvecserial -lsundials_sunlinsolklu -lsundials_sunmatrixsparse 


# Source files for binaries
APP_SRCS = $(SRC_DIR)/main.c $(SRC_DIR)/batt_model.c 


# Binary targets
BIN = $(BIN_DIR)/run_model


# All targets
.PHONY: all clean
all: $(BIN)


# Build binary
$(BIN): $(APP_SRCS)
	$(CC) $(CFLAGS) -o $(BIN) $(APP_SRCS) $(LIBS) 


# Compile source files to object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@


# Clean target
clean:
	rm -f $(BIN)
	rm -f dfn_gittr.csv 

