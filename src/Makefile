# CFLAGS=-Wall -g
# CFLAGS=-Wall -g -m64
# CFLAGS=-Wall -gdwarf-4 -g3
CFLAGS=-Wall
LINKER_FLAGS=-lm
CC=g++
TARGET=synder

OBJECTS=ia.o iv.o interval.o itree.o search.o arguments.o lev.o block.o contig.o genome.o contiguous_set.o synmap.o io.o analysis.o map.o main.o global.o

${TARGET}: ${OBJECTS}
	${CC} -o ${TARGET} ${OBJECTS} ${LINKER_FLAGS}

global.o: global.cpp global.hpp
	${CC} -c ${CFLAGS} global.cpp global.hpp

main.o: main.cpp
	${CC} -c ${CFLAGS} main.cpp

io.o: io.cpp io.hpp
	${CC} -c ${CFLAGS} io.cpp io.hpp

genome.o: genome.cpp genome.hpp
	${CC} -c ${CFLAGS} genome.cpp genome.hpp

analysis.o: analysis.cpp analysis.hpp
	${CC} -c ${CFLAGS} analysis.cpp analysis.hpp

contig.o: contig.cpp contig.hpp
	${CC} -c ${CFLAGS} contig.cpp contig.hpp

block.o: block.cpp block.hpp
	${CC} -c ${CFLAGS} block.cpp block.hpp

synmap.o: synmap.cpp synmap.hpp
	${CC} -c ${CFLAGS} synmap.cpp synmap.hpp

contiguous_set.o: contiguous_set.cpp contiguous_set.hpp
	${CC} -c ${CFLAGS} contiguous_set.cpp contiguous_set.hpp

arguments.o: arguments.cpp arguments.hpp version.hpp
	${CC} -c ${CFLAGS} arguments.cpp arguments.hpp

map.o: map.cpp map.hpp
	${CC} -c ${CFLAGS} map.cpp map.hpp

itree.o: itree/itree.cpp itree/itree.hpp
	${CC} -c ${CFLAGS} itree/itree.cpp itree/itree.hpp

interval.o: itree/interval.cpp itree/interval.hpp
	${CC} -c ${CFLAGS} itree/interval.cpp itree/interval.hpp

search.o: itree/search.cpp itree/search.hpp
	${CC} -c ${CFLAGS} itree/search.cpp itree/search.hpp

ia.o: itree/ia.cpp itree/ia.hpp
	${CC} -c ${CFLAGS} itree/ia.cpp itree/ia.hpp

iv.o: itree/iv.cpp itree/iv.hpp
	${CC} -c ${CFLAGS} itree/iv.cpp itree/iv.hpp

lev.o: lev.cpp lev.hpp
	${CC} -c ${CFLAGS} lev.cpp lev.hpp


profiler-build:
	${CC} -o ${TARGET} ${CFLAGS} -pg *.cpp *.hpp itree/*.cpp itree/*.hpp

clean:
	rm -f ${OBJECTS} *.{gch,o} itree/*.{gch,o} ${TARGET}
