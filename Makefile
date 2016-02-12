TARGET=synfull
LIBMAIN=src/main.a

${TARGET}: ${LIBMAIN}
	gcc -o ${TARGET} ${LIBMAIN}

${LIBMAIN}:
	make -C src --no-print-directory

docs:
	doxygen Doxyfile

clean:
	rm -f ${TARGET}

rclean:
	make clean
	make -C src --no-print-directory rclean
