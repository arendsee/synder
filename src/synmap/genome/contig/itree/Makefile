libitree.a: itree.o interval.o iv.o interval-tree.o ia.o
	ar cr libitree.a *.o

itree.o: itree.h itree.c iv.o interval-tree.o ia.o
	gcc -c itree.[hc]

interval-tree.o: interval-tree.c interval-tree.h ia.o
	gcc -c interval-tree.[hc]

iv.o: iv.h iv.c interval.o
	gcc -c iv.[hc]

ia.o: ia.c ia.h interval.o
	gcc -c ia.[hc]

interval.o: 
	gcc -c interval.[hc]

clean:
	rm -f  *.a *.o *.gch
