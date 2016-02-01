weave:
	cweave fagin.w && latex fagin.tex || rm fagin.tex
	rm -f fagin.{aux,dvi,idx,log,scn,toc} 

tangle:
	ctangle fagin.w && cc -o fagin fagin.c || rm fagin.c

clean:
	rm -f fagin.c fagin.tex fagin
