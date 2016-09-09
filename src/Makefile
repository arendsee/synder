# CFLAGS=-Wall -g
# CFLAGS=-Wall -g -m64
# CFLAGS=-Wall -gdwarf-4 -g3
CFLAGS=-Wall -gdwarf-4 -g3 -std=c++11
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
	lev.o              \
	block.o            \
	contig.o           \
	genome.o           \
	contiguous_set.o   \
	synmap.o           \
	io.o               \
	map.o              \
	main.o             \

${TARGET}: ${OBJECTS}
	${CC} -o ${TARGET} ${OBJECTS} ${LINKER_FLAGS}

%.o: %.cpp
	${CC} -c ${CFLAGS} $< -o $@
	
clean:
	rm -f ${OBJECTS} *.{gch,o} ${TARGET}
