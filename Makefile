# ============================================================
# Makefile for SMOLPACK — Smolyak Sparse-Grid Cubature Library
# ============================================================
#
# Build targets:
#   make          — build the test runner
#   make test     — build and run tests
#   make clean    — remove build artefacts
#
# Requirements: C11-capable compiler, math library (-lm)
# ============================================================

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -pedantic -O2
LDFLAGS  = -lm

SRC_DIR  = smolpack
TEST_DIR = tests
BUILD    = build

# ---- Source files ----
SRC_FILES = \
	$(SRC_DIR)/smolpack.c \
	$(SRC_DIR)/smolyak_init.c \
	$(SRC_DIR)/smolyak_weights.c \
	$(SRC_DIR)/smolyak_alg.c \
	$(SRC_DIR)/cc_init.c \
	$(SRC_DIR)/cc_weights.c \
	$(SRC_DIR)/cc_alg.c

TEST_FILES = \
	$(TEST_DIR)/genz.c \
	$(TEST_DIR)/test_smolpack.c

ALL_SRC = $(SRC_FILES) $(TEST_FILES)
OBJS    = $(patsubst %.c,$(BUILD)/%.o,$(ALL_SRC))
TARGET  = $(BUILD)/test_smolpack

# ---- Default target ----
all: $(TARGET)

$(TARGET): $(OBJS)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD)/$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -I$(SRC_DIR) -c -o $@ $<

$(BUILD)/$(TEST_DIR)/%.o: $(TEST_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -I$(SRC_DIR) -I$(TEST_DIR) -c -o $@ $<

# ---- Run tests ----
test: $(TARGET)
	./$(TARGET)

# ---- Clean ----
clean:
	rm -rf $(BUILD)

.PHONY: all test clean
