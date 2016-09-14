# # production build
CFLAGS=-Wall -O2 -std=c++11

# debug build
# CFLAGS=-Wall -gdwarf-4 -g3 -std=c++11

# # profiling build
# # NOTE: GNU profiler leaks memory, so turn off memory tests
# CFLAGS=-Wall -pg -std=c++11

LINKER_FLAGS=-lm
CC=g++
TARGET=synder

OBJECTS=               \
	global.o           \
	interval.o         \
	itree_result.o     \
	itree_node.o       \
	itree.o            \
	arguments.o        \
	block.o            \
	contig_result.o    \
	contig.o           \
	genome.o           \
	contiguous_set.o   \
	synmap.o           \
	map.o              \
	main.o             \

${TARGET}: ${OBJECTS}
	${CC} -o ${TARGET} ${OBJECTS} ${LINKER_FLAGS}

%.o: %.cpp
	${CC} -c ${CFLAGS} $< -o $@
	
clean:
	rm -f ${OBJECTS} *.{gch,o} ${TARGET}
