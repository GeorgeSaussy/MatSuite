matsuite: lib/mslib.c
	gcc -c -DBUILD_DLL lib/mslib.c

dox:
	cd lib; doxygen Doxyfile
	mv lib/html .
	mv lib/latex .

clean:
	$(RM) mslib.o
	$(RM) -rf html/
	$(RM) -rf latex/
