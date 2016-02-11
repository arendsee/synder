TARGET=libcontig.a
LIBITREE=itree/libitree.a
LIBBLOCK=block/libblock.a

${TARGET}: contig.o ${LIBBLOCK} ${LIBITREE}
	ar crTv ${TARGET} contig.o ${LIBBLOCK} ${LIBITREE}

contig.o: contig.h contig.h
	gcc -c contig.[ch]

${LIBBLOCK}:
	make -C block

${LIBITREE}:
	make -C itree

clean:
	rm -f ${TARGET} *.o *.gch

rclean:
	make clean
	make -C block clean
	make -C itree clean
