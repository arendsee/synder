TARGET=libcontig.a
LIBITREE=itree/libitree.a

.PHONY: itree-build ${TARGET}

itree-build:
	cd itree && ${MAKE}

${TARGET}: contig.o block.o ${LIBITREE}
	ar crT ${TARGET} contig.o block.o ${LIBITREE}

contig.o: contig.h contig.h
	gcc -c contig.[ch]

block.o:
	gcc -c block.[ch]

clean:
	rm -f ${TARGET} *.o *.gch

rclean:
	make clean
	make --no-print-directory -C itree clean
