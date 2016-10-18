matsuite: lib/mslib.c
	gcc -c -DBUILD_DLL lib/mslib.c

documentation:
	cd lib; doxygen Doxyfile

clean:
	$(RM) mslib.o
	$(RM) -rf lib/html/
	$(RM) -rf lib/latex/
    
