TARGET=synfull
LIBMAIN=src/main.a

${TARGET}: ${LIBMAIN}
	gcc -o ${TARGET} ${LIBMAIN}

${LIBMAIN}:
	make -C src

clean:
	rm -f ${TARGET}

rclean:
	make clean
	make -C src rclean
