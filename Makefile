TARGET=synder
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


.PHONY: uninstall
uninstall:
	rm -f ${PREFIX}/bin/synder


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
	./test/runtests.sh -xmdv -o log -a ${ARCHIVE} | tee ${ARCHIVE}/log

# Same as above, but does not test memory
.PHONY: dtest-leak
dtest-leak:
	${MAKE} ddclean
	mkdir ${ARCHIVE}
	./test/runtests.sh -xdv -o log -a ${ARCHIVE} | tee ${ARCHIVE}/log

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
	rm -f ${TEMP}

.PHONY: ddclean
ddclean:
	rm -f ${TEMP}
	rm -rf ${TEMP_DIR}
