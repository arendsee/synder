TARGET=libgenome.a
LIBCONTIG=contig/libcontig.a

${TARGET}: genome.o ${LIBCONTIG}
	ar crT ${TARGET} genome.o ${LIBCONTIG}

genome.o:
	gcc -c genome.[ch]

${LIBCONTIG}:
	make --no-print-directory -C contig

clean:
	rm -f ${TARGET} *.o *.gch

rclean:
	make clean
	make --no-print-directory -C contig rclean
