matsuite: lib/mslib.c
	gcc -c -DBUILD_DLL lib/mslib.c

dox:
	cd lib; doxygen Doxyfile
	mv lib/html .
	mv lib/latex .
	cd latex; make

tests:
	cp lib/mslib.c examples
	cp lib/mslib.h examples
	cd examples; g++ -o mult_test mslib.c testlib.cpp multtest.cpp -std=c++11
	$(RM) examples/mslib.c
	$(RM) examples/mslib.h

all:
	gcc -c -DBUILD_DLL lib/mslib.c
	cd lib; doxygen Doxyfile
	mv lib/html .
	mv lib/latex .
	cd latex; make
	cp lib/mslib.c examples
	cp lib/mslib.h examples
	cd examples; g++ -o mult_test mslib.c testlib.cpp multtest.cpp -std=c++11
	$(RM) examples/mslib.c
	$(RM) examples/mslib.h

clean:
	$(RM) mslib.o
	$(RM) -r html/
	$(RM) -r latex/
	$(RM) examples/mult_test
