TARGET=libsynmap.a
LIBGENOME=genome/libgenome.a
CFLAGS=-Wall -pedantic -g

${TARGET}: synmap.o ${LIBGENOME}
	ar crTv ${TARGET} synmap.o ${LIBGENOME}

synmap.o:
	gcc ${CFLAGS} -c synmap.[ch]

${LIBGENOME}:
	make -C genome

clean:
	rm -f ${TARGET} *.o *.gch

rclean:
	make clean
	make -C genome rclean
