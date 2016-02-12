libitree.a: itree.o interval.o iv.o interval-tree.o ia.o
	ar cr libitree.a *.o

itree.o: itree.h itree.c iv.h interval-tree.h ia.h
	gcc -c itree.[hc]

interval-tree.o: interval-tree.c interval-tree.h ia.h
	gcc -c interval-tree.[hc]

iv.o: iv.h iv.c interval.h
	gcc -c iv.[hc]

ia.o: ia.c ia.h interval.h
	gcc -c ia.[hc]

interval.o: interval.c interval.h
	gcc -c interval.[hc]

clean:
	rm -f  *.a *.o *.gch
