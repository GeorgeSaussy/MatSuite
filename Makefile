matsuite: lib/mslib.c
	gcc -c -DBUILD_DLL lib/mslib.c

dox:
	cd lib; doxygen Doxyfile
	mv lib/html .
	mv lib/latex .
	cd latex; make

clean:
	$(RM) mslib.o
	$(RM) -r html/
	$(RM) -r latex/
