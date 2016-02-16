CFLAGS=-Wall -Ofast
CC=gcc
TARGET=synfull

${TARGET}:
	${CC} -o ${TARGET} ${CFLAGS} *.c *.h itree/*.c itree/*.h

profiler-build:
	${CC} -o ${TARGET} ${CFLAGS} -pg *.c *.h itree/*.c itree/*.h

clean:
	rm -f *.o *.gch ${TARGET}
	rm -f itree/*.o itree/*.gch
