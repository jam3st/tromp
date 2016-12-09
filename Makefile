OPT   = -O3
FLAGS = -Wall -Wno-deprecated-declarations -D_POSIX_C_SOURCE=200112L $(OPT) -pthread 
GPP   = g++ -march=native -m64 -std=c++11 $(FLAGS)

all:	equi

equi:	equi.h equi_miner.h equi_miner.cpp Makefile
	$(GPP) -ggdb equi_miner.cpp blake/blake2b.cpp -o equi

clean:	
	-rm -f equi
