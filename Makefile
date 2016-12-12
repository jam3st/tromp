FLAGS = -Wall -Wno-deprecated-declarations  -march=native -m64 -std=c++11 -ggdb 

all:	equi

equi:	equi.h equi_miner.h equi_miner.cpp Makefile
	g++ -ggdb equi_miner.cpp blake/blake2b.cpp -o equi

clean:	
	-rm -f equi
