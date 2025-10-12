# Compiler settings
CC      = clang
AR      = ar
CFLAGS  = -Iinclude -I. -Wall -Wextra -Wpedantic -O3

# Directories
OBJ_DIR = obj
LIB     = geometry.a

# Local sources
SRC     = geometry.c

# External sources
EXT_INC = -Ipredicates/include
EXT_SRC = $(wildcard predicates/src/*.c)

# Combine all sources
ALL_SRC = $(SRC) $(EXT_SRC)

# Map sources to object files (flat obj directory)
ALL_OBJS = $(patsubst %.c,$(OBJ_DIR)/%.o,$(notdir $(ALL_SRC)))

# Default target
all: $(LIB)

# Pattern rule to build objects from sources
$(OBJ_DIR)/%.o: 
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) $(EXT_INC) -c $(filter %/$*.c $(notdir %/$*.c),$(ALL_SRC)) -o $@

# Build static library
$(LIB): $(ALL_OBJS)
	$(AR) rcs $@ $(ALL_OBJS)

# Clean
clean:
	rm -rf $(OBJ_DIR) $(LIB)

# Debug helper
print:
	@echo "Sources: $(ALL_SRC)"
	@echo "Objects: $(ALL_OBJS)"

