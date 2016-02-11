TARGET=libgenome.a
LIBCONTIG=contig/libcontig.a

${TARGET}: genome.o ${LIBCONTIG}
	ar crTv ${TARGET} genome.o ${LIBCONTIG}

genome.o:
	gcc -c genome.[ch]

${LIBCONTIG}:
	make -C contig

clean:
	rm -f ${TARGET} *.o *.gch

rclean:
	make clean
	make -C contig rclean
