TARGET=synfull
LIBMAIN=src/main.a

all: build-src ${TARGET}

build-src:
	cd src && ${MAKE}

${TARGET}:
	gcc -o ${TARGET} ${LIBMAIN}

docs:
	doxygen Doxyfile

clean:
	rm -f ${TARGET}

rclean:
	make clean
	make -C src --no-print-directory rclean
