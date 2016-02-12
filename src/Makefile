TARGET=main.a
LIBSYNMAP=synmap/libsynmap.a
OUI=ui/ui.o

.PHONY: ${LIBSYNMAP} ${TARGET} ${UIO}

${LIBSYNMAP}:
	cd synmap && ${MAKE}

${UIO}:
	cd ui && ${MAKE}

${TARGET}: main.o
	ar crT ${TARGET} main.o ${OUI} ${LIBSYNMAP}

main.o: main.c
	gcc -c main.c

clean:
	rm -f *.a *.o

rclean:
	make clean
	make --no-print-directory -C synmap rclean
	make --no-print-directory -C ui clean
