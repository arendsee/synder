TARGET=libcontig.a
LIBITREE=itree/libitree.a
LIBBLOCK=block/libblock.a

${TARGET}: contig.o ${LIBBLOCK} ${LIBITREE}
	ar crT ${TARGET} contig.o ${LIBBLOCK} ${LIBITREE}

contig.o: contig.h contig.h
	gcc -c contig.[ch]

${LIBBLOCK}:
	make --no-print-directory -C block

${LIBITREE}:
	make --no-print-directory -C itree

clean:
	rm -f ${TARGET} *.o *.gch

rclean:
	make clean
	make --no-print-directory -C block clean
	make --no-print-directory -C itree clean
