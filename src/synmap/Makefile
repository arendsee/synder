TARGET=libsynmap.a
LIBGENOME=genome/libgenome.a
CFLAGS=-Wall -pedantic -g

${TARGET}: synmap.o synmap_io.o synmap_analysis.o ${LIBGENOME}
	ar crT ${TARGET} synmap.o synmap_io.o synmap_analysis.o ${LIBGENOME}

synmap.o:
	gcc ${CFLAGS} -c synmap.[ch]

synmap_io.o:
	gcc ${CFLAGS} -c synmap_io.[ch]
	
synmap_analysis.o:
	gcc ${CFLAGS} -c synmap_analysis.[ch]

${LIBGENOME}:
	make -C genome

clean:
	rm -f ${TARGET} *.o *.gch

rclean:
	make clean
	make -C genome rclean
