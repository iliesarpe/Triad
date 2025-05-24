CC = g++
CFLAGS = -v --std=c++2a -g -O3 -fopenmp -march=native#-fsanitize=address 
GUROBI_INCL="" #FIXME: ADD YOUT PATH
GUROBI_LIB="" #FIXME: ADD YOUT PATH
FOLLIBS = ${PWD}/libs/
EXEC = clustcoeff


all: main.cpp ${FOLLIBS}structures.h structures.o
	$(CC) ${CFLAGS} main.cpp structures.o -o ${EXEC} -Ilibs/ -I./statLibs/ -I${GUROBI_INCL} -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -L${GUROBI_LIB} -lgurobi_c++ -lhdf5 -ldl -lgurobi110

structures.o: ${FOLLIBS}structures.cpp
	$(CC) -c ${CFLAGS} ${FOLLIBS}structures.cpp -Ilibs/ -I./statLibs/ -I${GUROBI_INCL} -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -L${GUROBI_LIB}  -lgurobi_c++ -lhdf5 -ldl -lgurobi110

clean:
	rm ${EXEC} *.o
	rm -rf *.dSYM
