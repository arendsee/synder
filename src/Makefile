CFLAGS=-Wall -g -O2
CC=gcc
TARGET=run
OBJECTS=main.o io.o genome.o analysis.o contig.o itree/interval.o itree/itree.o itree/search.o itree/ia.o itree/iv.o block.o synmap.o ui.o

${TARGET}: ${OBJECTS}
	${CC} -o ${TARGET} ${CFLAGS} ${OBJECTS}



profiler-build:
	${CC} -o ${TARGET} ${CFLAGS} -pg *.c *.h itree/*.c itree/*.h

clean:
	rm -f ${OBJECTS} *.gch itree/*.gch ${TARGET}
