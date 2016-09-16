# # production build
# CFLAGS=-Wall -O2 -std=c++11

# debug build
CFLAGS=-Wall -pedantic -gdwarf-4 -g3 -std=c++11

# # profiling build
# # NOTE: GNU profiler leaks memory, so turn off memory tests
# CFLAGS=-Wall -pg -std=c++11

LINKER_FLAGS=-lm
CC=g++
TARGET=synder

OBJECTS=               \
	global.o           \
	interval.o         \
	bound.o            \
	interval_tree.o    \
	arguments.o        \
	block.o            \
	contig.o           \
	genome.o           \
	contiguous_set.o   \
	synmap.o           \
	main.o             \

${TARGET}: ${OBJECTS}
	${CC} -o ${TARGET} ${OBJECTS} ${LINKER_FLAGS}

%.o: %.cpp
	${CC} -c ${CFLAGS} $< -o $@
	
clean:
	rm -f ${OBJECTS} *.{gch,o} ${TARGET}
