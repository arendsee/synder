CFLAGS=-Wall -g -O2
CC=gcc
TARGET=run
OBJECTS=main.o test.o io.o genome.o analysis.o contig.o interval.o itree.o search.o ia.o iv.o block.o synmap.o ui.o contiguous.o

${TARGET}: ${OBJECTS}
	${CC} -o ${TARGET} ${CFLAGS} ${OBJECTS}

main.o: main.c
	${CC} -c ${CFLAGS} main.c

test.o: test.c test.h
	${CC} -c ${CFLAGS} test.[ch]

io.o: io.c io.h
	${CC} -c ${CFLAGS} io.[ch]

genome.o: genome.c genome.h
	${CC} -c ${CFLAGS} genome.[ch]

analysis.o: analysis.c analysis.h
	${CC} -c ${CFLAGS} analysis.[ch]

contig.o: contig.c contig.h
	${CC} -c ${CFLAGS} contig.[ch]

block.o: block.c block.h
	${CC} -c ${CFLAGS} block.[ch]

synmap.o: synmap.c synmap.h
	${CC} -c ${CFLAGS} synmap.[ch]

ui.o: ui.c ui.h
	${CC} -c ${CFLAGS} ui.[ch]

contiguous.o: contiguous.c contiguous.h
	${CC} -c ${CFLAGS} contiguous.[ch]

itree.o: itree/itree.c itree/itree.h
	${CC} -c ${CFLAGS} itree/itree.[ch]

interval.o: itree/interval.c itree/interval.h
	${CC} -c ${CFLAGS} itree/interval.[ch]

search.o: itree/search.c itree/search.h
	${CC} -c ${CFLAGS} itree/search.[ch]

ia.o: itree/ia.c itree/ia.h
	${CC} -c ${CFLAGS} itree/ia.[ch]

iv.o: itree/iv.c itree/iv.h
	${CC} -c ${CFLAGS} itree/iv.[ch]



profiler-build:
	${CC} -o ${TARGET} ${CFLAGS} -pg *.c *.h itree/*.c itree/*.h

clean:
	rm -f ${OBJECTS} *.{gch,o} itree/*.{gch,o} ${TARGET}
