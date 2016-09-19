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

OBJECTS=                   \
	global.o               \
	arguments.o            \
	interval.o             \
	linked_interval.o      \
	bound.o                \
	interval_tree.o        \
	block.o                \
	contiguous_set.o       \
	interval_set.o         \
	many_blocks.o          \
	contig.o               \
	genome.o               \
	synmap.o               \
	main.o                 \

${TARGET}: ${OBJECTS}
	${CC} -o ${TARGET} ${OBJECTS} ${LINKER_FLAGS}

%.o: %.cpp
	${CC} -c ${CFLAGS} $< -o $@
	
clean:
	rm -f ${OBJECTS} *.{gch,o} ${TARGET}
