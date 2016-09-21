# # production build
# CFLAGS=-Wall -O2 -std=c++11

# debug build
CFLAGS=-Wall -pedantic -gdwarf-4 -g3 -std=c++11
INSUFFERABLE=-Wshadow # -Wextra -Wold-style-cast -Wcast-align -Weffc++ -Wmisleading-indentation

# # profiling build
# # NOTE: GNU profiler leaks memory, so turn off memory tests
# CFLAGS=-Wall -pg -std=c++11

LINKER_FLAGS=-lm
CC=g++
TARGET=synder

OBJECTS=                   \
	global.o               \
	arguments.o            \
	bound.o                \
	block.o                \
	contiguous_set.o       \
	many_blocks.o          \
	many_contiguous_sets.o \
	contig.o               \
	genome.o               \
	synmap.o               \
	main.o                 \

${TARGET}: ${OBJECTS}
	${CC} -o ${TARGET} ${OBJECTS} ${LINKER_FLAGS}

%.o: %.cpp
	${CC} -c ${CFLAGS} ${INSUFFERABLE} $< -o $@
	
clean:
	rm -f ${OBJECTS} *.{gch,o} ${TARGET}
