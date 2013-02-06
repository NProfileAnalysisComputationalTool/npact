FLAGS = -ggdb
ENDFLAGS=  #it doesn't work unless this is a the end of the command line.
SRCDIR = src
OBJDIR = obj
BINDIR = pynpact/bin
CC = gcc

SRC = $(shell find $(SRCDIR) | grep .c)
OBJ = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRC))
_PROGS = acgt_gamma AG Allplots atg extract gtg nprofile stops ttg
PROGS = $(patsubst %,$(BINDIR)/%,$(_PROGS))

DEPS = $(shell find $(SRCDIR) | grep .h)
DEPSO = $(pathsubs %,$(OBJDIR)/%.o,$(_DEPS))

.PHONY: all clean init

all: init $(PROGS)
clean:
	rm -rf $(OBJDIR) $(BINDIR)

init:
	mkdir -p $(OBJDIR) $(BINDIR)

#A general rule of make an object file for any C file, rebuild if the
#corresponding C file changes or any of the include files change
#(heavy handed but easy).
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) $(FLAGS) -c $< -o $@


#Link acgt_gamma if the object file changes include lib m
$(BINDIR)/acgt_gamma: $(OBJDIR)/acgt_gamma.o $(OBJDIR)/util.o
	$(CC) $(FLAGS) -o $@ $^  -lm

#Allplots needs util and lib m
$(BINDIR)/Allplots: $(OBJDIR)/Allplots.o $(OBJDIR)/util.o
	$(CC) $(FLAGS) -o $@ $^  -lm

#anything else can be linked with just its own object file.
$(BINDIR)/%: $(OBJDIR)/%.o
	$(CC) $(FLAGS) -o $@ $^  $(ENDFLAGS)

