# project Makefile
# Use gcc compiler
CC := gcc

# OS select
ifeq ($(OS),Windows_NT)
	OBJS2 = $(subst /,\,$(OBJS))
    RM = del /Q /F
endif

# directory path
LIBDIR = ./lib
OBJDIR = ../obj
INC = -I../inc

# sorce code
TARGET := ../bin/exe
MAIN = ILU0_bicgstab.c
LIBS = $(wildcard $(LIBDIR)/*.c)
SRCS = $(MAIN) $(LIBS)
OBJS = $(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(SRCS)))

# compiler option
CFLAGS := -c
LDFLAGS := 
RM = rm -f

$(TARGET) : $(OBJS)
	$(CC) $^ -o $@ $(LDFLAGS)

$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) $(INC) $< -o $@ 

# all
.PHONY: all
all: $(TARGET)

# debug
.PHONY: debug
debug: CFLAGS += -g -O0
debug: all

# release
.PHONY: release
release: CFLAGS += -o2 -s
release: all

# parallel
.PHONY: parallel
parallel: CC = icc
parallel: CFLAGS += -qopenmp
parallel: LDFLAGS += -parallel
parallel: all

# make clean
.PHONY: clean
clean:
	$(RM) $(OBJS)

# make run
.PHONY: run
run:
	../bin/exe

