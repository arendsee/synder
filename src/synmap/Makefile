TARGET=libsynmap.a
LIBCONTIG=contig/libcontig.a
CFLAGS=-Wall -pedantic -g

.PHONY: ${LIBCONTIG} ${TARGET}

${LIBCONTIG}:
	cd contig && ${MAKE}

${TARGET}: synmap.o synmap_io.o synmap_analysis.o genome.o
	ar crT ${TARGET} synmap.o synmap_io.o synmap_analysis.o ${LIBCONTIG}

synmap.o:
	gcc ${CFLAGS} -c synmap.[ch]

synmap_io.o:
	gcc ${CFLAGS} -c synmap_io.[ch]
	
synmap_analysis.o:
	gcc ${CFLAGS} -c synmap_analysis.[ch]

genome.o:
	gcc ${CFLAGS} -c genome.[ch]

clean:
	rm -f ${TARGET} *.o *.gch

rclean:
	make clean
	make -C contig rclean
