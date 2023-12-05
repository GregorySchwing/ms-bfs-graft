CC = g++ #-std=c99
cc = gcc
CFLAGS = -fopenmp -Drestrict=__restrict__ -O2 -DNDEBUG -ffast-math # -g -pg
LDFLAGS = -O2

all: msBFSGraft msBFSGraft_lib.a

graphgenBP.o: graphgenBP.h graphgenBP.cpp ThreadedMMReader.h
	$(CC) $(CFLAGS) -c -o graphgenBP.o graphgenBP.cpp 

msBFSGraft.o: msBFSGraft.cpp graphgenBP.h graphgenBP.o maximalMatching.o 
	$(CC) $(CFLAGS) -c -o msBFSGraft.o msBFSGraft.cpp 


pf.o: PothenFan.cpp graphgenBP.h graphgenBP.o maximalMatching.o
	$(CC) $(CFLAGS) -c -o pf.o PothenFan.cpp 

mmio.o: mmio.c
	$(CC) $(CFLAGS) -Wno-write-strings -c -o mmio.o mmio.c

msBFSGraft: msBFSGraft.o maximalMatching.o mmio.o 
	$(CC) $(CFLAGS) $(LDFLAGS)  -o msBFSGraft msBFSGraft.o graphgenBP.o mmio.o maximalMatching.o -lm -lstdc++

pf: pf.o  maximalMatching.o
	$(CC) $(CFLAGS) $(LDFLAGS)  -o pf pf.o graphgenBP.o maximalMatching.o -lm -lstdc++

maximalMatching.o: maximalMatching.cpp maximalMatching.h
	$(CC) $(CFLAGS) -c -o maximalMatching.o maximalMatching.cpp 


msBFSGraft_lib.o: msBFSGraft_lib.cpp graphgenBP.h graphgenBP.o maximalMatching.o 
	$(CC) $(CFLAGS) -c -o msBFSGraft_lib.o msBFSGraft_lib.cpp 

msBFSGraft_lib.a: msBFSGraft_lib.o maximalMatching.o mmio.o 
	ar rvs  msBFSGraft_lib.a msBFSGraft_lib.o graphgenBP.o mmio.o maximalMatching.o

clean:
	-rm -f msBFSGraft *.o
