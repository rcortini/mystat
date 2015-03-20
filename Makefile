program=mystat
program_sources=
program_main=$(program).c
program_objects=$(program_sources:.c=.o)
program_dependencies=

# flags
CFLAGS=-Wall -Wextra -Wshadow -g --pedantic -O0 $(shell pkg-config mylib --cflags)
LDFLAGS=$(shell pkg-config mylib --libs)

# programs
RM_F=rm -f
TAR=tar
CC=gcc

# targets
.PHONY: clean tarball install

all : $(program)

$(program): $(program_objects) $(program_main)
	$(CC) $(CFLAGS) -o $(program) $(program_main) $(program_objects) $(LDFLAGS)

# object files
%.o: %.c $(dependencies)
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

# clean targets
clean :
	$(RM_F) $(program_objects) $(program)

# install
install : $(program)
	cp $(program) ${SOFT_DIR}/bin

# tarball
tarball :
	$(TAR) -cf $(program).tar $(program_sources) $(program_main) $(program_dependencies) Makefile
