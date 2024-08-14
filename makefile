#!MAKE
SRCDIR = src
OBJDIR = obj
CC = gcc
CFLAGS = -fopenmp -g -lc
ASAN_FLAGS = -fsanitize=address -fno-omit-frame-pointer -Wno-format-security -static-libasan

OBJS = $(OBJDIR)/mesh.o
OBJS += $(OBJDIR)/matrix.o
OBJS += $(OBJDIR)/params.o
OBJS += $(OBJDIR)/assemble.o

diffusion: $(OBJS)
	$(CC) -o $@ $(SRCDIR)/main.c $^ $(CFLAGS) $(ASAN_FLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(SRCDIR)/%.h
	$(CC) -c -o $@ $< $(CFLAGS) $(ASAN_FLAGS)

clean:
	rm diffusion $(OBJDIR)/*.o
