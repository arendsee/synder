TARGET=main.a
LIBSYNMAP=synmap/libsynmap.a
OUI=ui/ui.o

${TARGET}: main.o ${OUI} ${LIBSYNMAP}
	ar crT ${TARGET} main.o ${OUI} ${LIBSYNMAP}

main.o: main.c
	gcc -c main.c

${LIBSYNMAP}:
	make --no-print-directory -C synmap

${UIO}:
	make --no-print-directory -C ui

clean:
	rm -f *.a *.o

rclean:
	make clean
	make --no-print-directory -C synmap rclean
	make --no-print-directory -C ui clean
