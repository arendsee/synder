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


.PHONY: docs
docs:
	doxygen Doxyfile


# ===================================================================
# Convenience functions

.PHONY: clean
clean:
	rm -f ${TARGET}
	rm -f vgcore.* gmon.out *log tags valgrind*
	rm -rf zzz* db d e g x c o m v z gdb.txt
	cd src && ${MAKE} clean


# Debug test
# * stops at first error
# * preps for GDB-based debugging
# * links valgrind output
.PHONY: dtest
dtest:
	./test/runtests.sh -xdo `tty`
	# Compile ctags, no worries if it fails
	ctags . 2> /dev/null

.PHONY: dclean
dclean:
	rm -f vgcore.* gmon.out *log valgrind*
	rm -rf zzz* db d e g x c o m v z gdb.txt


# Runs the sample data, linking files for review
.PHONY: sample
sample:
	${MAKE}
	./synder -d sample-inputs/a-b.syn a b db
	ln -sf sample-inputs/a.gff g
	ln -sf db/a_b.txt d
	./synder -i g -s d -c search
