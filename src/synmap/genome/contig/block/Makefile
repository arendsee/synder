TARGET=libblock.a

${TARGET}: block.o
	ar cr ${TARGET} block.o

block.o: block.c block.h
	cc -c block.c block.h

clean:
	rm -f ${TARGET} *.o *.gch
