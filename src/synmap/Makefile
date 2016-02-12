TARGET=libsynmap.a
LIBCONTIG=contig/libcontig.a
CFLAGS=-Wall -pedantic -g
OBJECTS=synmap.o synmap_io.o synmap_analysis.o genome.o

all: make-contig ${TARGET}

make-contig:
	cd contig && ${MAKE}

${TARGET}: ${OBJECTS}
	ar crT ${TARGET} ${OBJECTS} ${LIBCONTIG}

synmap.o: synmap.c synmap.h genome.h
	gcc ${CFLAGS} -c synmap.[ch]

synmap_io.o: synmap_io.c synmap_io.h genome.h
	gcc ${CFLAGS} -c synmap_io.[ch]
	
synmap_analysis.o: synmap_analysis.c synmap_analysis.h genome.h
	gcc ${CFLAGS} -c synmap_analysis.[ch]

genome.o: genome.c genome.h
	gcc ${CFLAGS} -c genome.[ch]

clean:
	rm -f ${TARGET} ${OBJECTS} *.gch

rclean:
	make clean
	make -C contig rclean
