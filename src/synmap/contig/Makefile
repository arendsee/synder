TARGET=libcontig.a
LIBITREE=itree/libitree.a

all: make-itree ${TARGET}

make-itree:
	cd itree && ${MAKE}

${TARGET}: contig.o block.o
	ar crT ${TARGET} contig.o block.o ${LIBITREE}

contig.o: contig.c contig.h block.h
	gcc -c contig.[ch]

block.o: block.c block.h
	gcc -c block.[ch]

clean:
	rm -f ${TARGET} *.o *.gch

rclean:
	make clean
	make --no-print-directory -C itree clean
