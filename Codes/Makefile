# Makefile for hexatic_g6_avg pipeline
# Usage:
#   make            # release build
#   make DEBUG=1    # debug build (-g, -O0)
#   make clean
#   make run ARGS="..."   # run the program with ARGS
#
# Triangle handling:
# - If triangle.c exists in the project (root or src/), it will be compiled automatically.
# - Otherwise set TRIANGLE_LIB to a library (e.g. -ltriangle) and TRIANGLE_INC if header is elsewhere.

# ---------- Config ----------
CC       ?= gcc
CFLAGS   ?= -std=c99 -O3 -Wall -Wextra -pipe
DEBUG_CFLAGS = -g -O0 -DDEBUG -fsanitize=address,undefined

# If you have a triangle library installed, set TRIANGLE_LIB (e.g. -ltriangle)
TRIANGLE_LIB ?=
# If triangle.h is in a non-standard include dir, set this (e.g. -I/usr/local/include)
TRIANGLE_INC ?=

# Source layout (edit if your files live elsewhere)
SRCDIR   := .
SRCS     := $(SRCDIR)/main.c \
            $(SRCDIR)/utils.c \
            $(SRCDIR)/io.c \
            $(SRCDIR)/clusters.c \
            $(SRCDIR)/com.c \
            $(SRCDIR)/delaunay.c \
            $(SRCDIR)/psi6.c \
            $(SRCDIR)/g6accum.c

# If triangle.c is present in project, compile it
TRI_CANDIDATES := triangle.c 
TRI_SRC := $(firstword $(filter %,$(wildcard $(TRI_CANDIDATES))))
ifeq ($(TRI_SRC),)
# No triangle.c found in project. Will try to link using TRIANGLE_LIB (user must set)
TRI_OBJ :=
else
TRI_OBJ := $(TRI_SRC:.c=.o)
SRCS   += $(TRI_SRC)
endif

# Output program
PROG := hexatic_g6_avg

# Derived
OBJS := $(SRCS:.c=.o)
DEPS := $(OBJS:.o=.d)

# Allow overriding compiler flags (e.g. add -I)
# CPPFLAGS ?= $(TRIANGLE_INC)
# CPPFLAGS ?= $(TRIANGLE_INC) -DTRILIBRARY
# AFTER
CPPFLAGS ?= $(TRIANGLE_INC) -DTRILIBRARY -DNOTIMING

# If DEBUG, adjust flags
ifeq ($(DEBUG),1)
CFLAGS := $(CFLAGS) $(DEBUG_CFLAGS)
endif

LDFLAGS ?=
LDLIBS  := $(TRIANGLE_LIB) -lm

.PHONY: all clean run

all: $(PROG)

# Link
$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Compile C -> object with dependency generation
# -MMD -MP creates .d files for header deps
%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -MMD -MP -c -o $@ $<

# Include generated dependency files
-include $(DEPS)

# Convenience run target: make run ARGS="..." 
run: $(PROG)
	./$(PROG) $(ARGS)

# Clean
clean:
	$(RM) $(PROG) $(OBJS) $(DEPS)

# Show configuration
info:
	@echo "CC       = $(CC)"
	@echo "CFLAGS   = $(CFLAGS)"
	@echo "CPPFLAGS = $(CPPFLAGS)"
	@echo "LDFLAGS  = $(LDFLAGS)"
	@echo "LDLIBS   = $(LDLIBS)"
	@echo "SRCS     = $(SRCS)"
	@echo "OBJS     = $(OBJS)"
	@echo "TRI_SRC  = $(TRI_SRC)"
