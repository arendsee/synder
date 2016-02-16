TARGET=synfull

all:
	cd src && ${MAKE}
	mv src/${TARGET} ${TARGET}

docs:
	doxygen Doxyfile

clean:
	rm -f ${TARGET}

rclean:
	make clean
	make -C src --no-print-directory clean
