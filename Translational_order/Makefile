CC       ?= gcc
CFLAGS   ?= -std=c99 -O3 -Wall -Wextra -pipe
DEBUG_CFLAGS = -g -O0 -DDEBUG -fsanitize=address,undefined

TRIANGLE_LIB ?=
TRIANGLE_INC ?=

SRCS := main.c utils.c io.c clusters.c com.c delaunay.c triangle.c psi6.c graccum.c gtaccum.c
OBJS := $(SRCS:.c=.o)
DEPS := $(OBJS:.o=.d)
PROG := translational_order

CPPFLAGS ?= $(TRIANGLE_INC) -DTRILIBRARY -DNOTIMING
LDFLAGS ?=
LDLIBS := $(TRIANGLE_LIB) -lm

ifeq ($(DEBUG),1)
CFLAGS := $(CFLAGS) $(DEBUG_CFLAGS)
endif

.PHONY: all clean
all: $(PROG)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -MMD -MP -c -o $@ $<

-include $(DEPS)

clean:
	$(RM) $(PROG) $(OBJS) $(DEPS)
