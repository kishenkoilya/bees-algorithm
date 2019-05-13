all:
	c++ -I../ -I../utils tutorial.cpp -o tutorial -fopenmp
clean:
	rm -rf *.o tutorial ht
