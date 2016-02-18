TARGET=run

${TARGET}: ${OBJECTS}
	cd src && ${MAKE}
	mv src/${TARGET} ${TARGET}

docs:
	doxygen Doxyfile

clean:
	rm -f ${TARGET}
	cd src && ${MAKE} clean
