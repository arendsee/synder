TARGET=synder
DBSCRIPT=make-synder-db.sh
DBSCRIPT_DIR=util
PREFIX=/usr/local

all:
	cd src && ${MAKE}
	mv src/${TARGET} ${TARGET}


.PHONY: test 
test:
	./test/runtests.sh


.PHONY: install
install:
	mkdir -p ${PREFIX}/bin
	cp -f ${TARGET} ${PREFIX}/bin
	cp -f ${DBSCRIPT_DIR}/${DBSCRIPT} ${PREFIX}/bin


.PHONY: uninstall
uninstall:
	rm -f ${PREFIX}/bin/synder
	rm -f ${PREFIX}/bin/${DBSCRIPT}


# ===================================================================
# Documentation

.PHONY: docs
docs:
	doxygen Doxyfile


# ===================================================================
# Utilities

.PHONY: tags
tag:
	ctags .


# ===================================================================
# Debugging
ARCHIVE=ark

# Debug test
# * stops at first error
# * preps for GDB-based debugging
# * links valgrind output
.PHONY: dtest
dtest:
	${MAKE} ddclean
	mkdir ${ARCHIVE}
	./test/runtests.sh -mdv -o log -a ${ARCHIVE} | tee ${ARCHIVE}/log

# Same as above, but does not test memory
.PHONY: dtest-leak
dtest-leak:
	${MAKE} ddclean
	mkdir ${ARCHIVE}
	./test/runtests.sh -dv -o log -a ${ARCHIVE} | tee ${ARCHIVE}/log

# Runs the sample data, linking files for review
.PHONY: sample
sample:
	${MAKE} ddclean
	${MAKE}
	./synder -d sample-inputs/a-b.syn a b db
	ln -sf sample-inputs/a.gff g
	ln -sf db/a_b.txt d
	./synder -i g -s d -c search


# ===================================================================
# Cleaning functions

TEMP=vgcor* gmon.out *log valgrind* gdb.txt [a-z] zz* expected-output gdb input.gff observed-output synteny-map.tab synmap.txt .gdb_cmd run
TEMP_DIR=db doc ${ARCHIVE}

.PHONY: clean
clean:
	rm -f ${TARGET} ${TEMP}
	rm -rf ${TEMP_DIR}
	cd src && ${MAKE} clean

.PHONY: dclean
dclean:
	rm -rf db
	rm -f ${TEMP}

.PHONY: ddclean
ddclean:
	rm -f ${TEMP}
	rm -rf ${TEMP_DIR}
