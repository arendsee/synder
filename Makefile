TARGET=synder
DBSCRIPT=make-synder-db.sh
DBSCRIPT_DIR=util
PREFIX=/usr/local

all:
	cd src && ${MAKE}
	mv src/${TARGET} ${TARGET}

docs:
	doxygen Doxyfile

clean:
	rm -f ${TARGET}
	cd src && ${MAKE} clean

install:
	mkdir -p ${PREFIX}/bin
	cp -f ${TARGET} ${PREFIX}/bin
	cp -f ${DBSCRIPT_DIR}/${DBSCRIPT} ${PREFIX}/bin

uninstall:
	rm -f ${PREFIX}/bin/synder
	rm -f ${PREFIX}/bin/${DBSCRIPT}
