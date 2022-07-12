# project Makefile
# Use gcc compiler
CC := gcc

# directory path
LIBDIR = lib
OBJDIR = ../obj
INC = -I../inc

# sorce
TARGET := ../bin/exe
MAIN = ILU0_bicgstab.c
LIBS = $(wildcard $(LIBDIR)/*.c)
SRCS = $(MAIN) $(LIBS)
OBJS = $(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(SRCS)))
OBJSW = $(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(SRCS)))

# compiler option
CFLAGS := -c
LDFLAGS := 
RM = rm -f

# OS select
ifeq ($(OS),Windows_NT)
	OBJSW = $(subst /,\,$(OBJS))
    RM = del /Q /F
endif

$(TARGET) : $(OBJS)
	$(CC) $^ -o $@ $(LDFLAGS)

$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) $(INC) $< -o $@ 

.PHONY: all debug release parallel clean run printenv

# all
all: $(TARGET)

# debug
debug: CFLAGS += -g -O0
debug: all

# release
release: CFLAGS += -o2 -s
release: all

# parallel
parallel: CC = icc
parallel: CFLAGS += -qopenmp
parallel: LDFLAGS += -parallel
parallel: all

# make clean
clean:
	$(RM) $(OBJSW)

# make run
run:
	../bin/exe

printenv:
	@echo OS         : $(OS)
	@echo CC         : $(CC)
	@echo OBJS       : $(OBJS)
	@echo OBJSW      : $(OBJSW)
	@echo RM         : $(RM)
	@echo LODIR      : $(LODIR)
	@echo REDIR      : $(REDIR)
	