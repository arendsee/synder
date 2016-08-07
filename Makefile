TARGET=synder
DBSCRIPT=make-synder-db.sh
DBSCRIPT_DIR=util
PREFIX=/usr/local

all:
	cd src && ${MAKE}
	mv src/${TARGET} ${TARGET}

docs:
	doxygen Doxyfile

install:
	mkdir -p ${PREFIX}/bin
	cp -f ${TARGET} ${PREFIX}/bin
	cp -f ${DBSCRIPT_DIR}/${DBSCRIPT} ${PREFIX}/bin

uninstall:
	rm -f ${PREFIX}/bin/synder
	rm -f ${PREFIX}/bin/${DBSCRIPT}

.PHONY: clean
clean:
	rm -f ${TARGET}
	rm -f vgcore.*
	cd src && ${MAKE} clean

.PHONY: test 
test:
	./test/runtests.sh
