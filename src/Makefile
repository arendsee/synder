TARGET=main.a
LIBSYNMAP=synmap/libsynmap.a
OUI=ui/ui.o

${TARGET}: main.o ${OUI} ${LIBSYNMAP}
	ar crTv ${TARGET} main.o ${OUI} ${LIBSYNMAP}

main.o: main.c
	gcc -c main.c

${LIBSYNMAP}:
	make -C synmap

${UIO}:
	make -C ui

clean:
	rm -f *.a *.o

rclean:
	make clean
	make -C synmap rclean
	make -C ui clean
