FLAGS = -Wall -Wno-deprecated-declarations  -march=native -m64 -std=c++14  -ggdb # -O3 #  -ggdb
all:	equi

equi:	equi.h equi_miner.h equi_miner.cpp blake2b.h basetype.h equihash.h Makefile
	g++ ${FLAGS} equi_miner.cpp blake/blake2b.cpp -o equi

clean:	
	-rm -f equi
