TARGET=synder
DBSCRIPT=make-synder-db.sh
DBSCRIPT_DIR=util
PREFIX=/usr/local

all:
	cd src && ${MAKE}
	mv src/${TARGET} ${TARGET}

.PHONY: sample
sample:
	${MAKE}
	./synder -d sample-inputs/a-b.syn a b db
	ln -sf sample-inputs/a.gff g
	ln -sf db/a_b.txt d
	./synder -i g -s d -c search

.PHONY: clean-sample
clean-sample:
	rm -rf g d db

.PHONY: docs
docs:
	doxygen Doxyfile

.PHONY: install
install:
	mkdir -p ${PREFIX}/bin
	cp -f ${TARGET} ${PREFIX}/bin
	cp -f ${DBSCRIPT_DIR}/${DBSCRIPT} ${PREFIX}/bin

.PHONY: uninstall
uninstall:
	rm -f ${PREFIX}/bin/synder
	rm -f ${PREFIX}/bin/${DBSCRIPT}

.PHONY: clean
clean:
	rm -f ${TARGET}
	rm -f vgcore.* gmon.out *log tags valgrind*
	rm -rf zzz* db d e g
	cd src && ${MAKE} clean

.PHONY: test 
test:
	./test/runtests.sh
