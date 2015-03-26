CC = gcc
CFLAGS = -g -c -Wall
OBJS = build/profiles.o build/paramprof.o build/bheap.o build/idr.o build/alignio.o build/xcorr.o build/iofile.o build/paramclust.o build/cluster.o build/hierarchical.o build/itvltree.o build/dtw.o build/profilemap.o build/annotate.o
TESTOBJS = build/tcparamprof.o

all : srnap


# Build and compile executables
 
srnap : profiles.o annotate.o srnap.o
	$(CC) -o bin/srnap $(OBJS) build/srnap.o -Llib/ -lgsl -lgslcblas -lm -lz -lpthread -lbam 

srnap.o : setup
	$(CC) $(CFLAGS) src/srnap.c -Isrc/include -o build/srnap.o

# Unit test

utest: profiles.o annotate.o utest.o
	gcc -o bin/utest -L/home/apages/tools/lcut-0.3.0/lib/ -Llib/ $(OBJS) $(TESTOBJS) build/utest.o -llcut -lgsl -lgslcblas -lm -lz -lpthread -lbam

utest.o: tcparamprof.o
	$(CC) $(CFLAGS) test/src/utest/utest.c -Isrc/include/ -Itest/src/utest/include -I/home/apages/tools/lcut-0.3.0/include/ -o build/utest.o


# Compile shared objects

annotate.o : paramclust.o xcorr.o iofile.o dtw.o hierarchical.o profilemap.o
	$(CC) $(CFLAGS) src/annotate/annotate.c -Isrc/include -o build/annotate.o

profilemap.o : itvltree.o
	$(CC) $(CFLAGS) src/annotate/profilemap.c -Isrc/include/ -o build/profilemap.o

hierarchical.o : cluster.o
	$(CC) $(CFLAGS) src/annotate/hierarchical.c -Isrc/include -o build/hierarchical.o

dtw.o : setup
	$(CC) $(CFLAGS) src/annotate/dtw.c -Isrc/include -o build/dtw.o

cluster.o : setup
	$(CC) $(CFLAGS) src/annotate/cluster.c -Isrc/include -o build/cluster.o

paramclust.o : setup
	$(CC) $(CFLAGS) src/annotate/paramclust.c -Isrc/include -o build/paramclust.o

xcorr.o : setup
	$(CC) $(CFLAGS) src/annotate/xcorr.c -Isrc/include/ -o build/xcorr.o

iofile.o : setup
	$(CC) $(CFLAGS) src/annotate/iofile.c -Isrc/include/ -o build/iofile.o

itvltree.o : setup
	$(CC) $(CFLAGS) src/annotate/itvltree.c -Isrc/include/ -o build/itvltree.o

profiles.o : paramprof.o bheap.o idr.o alignio.o
	$(CC) $(CFLAGS) src/profiles/profiles.c -Isrc/include -o build/profiles.o

paramprof.o : setup
	$(CC) $(CFLAGS) src/profiles/paramprof.c -Isrc/include -o build/paramprof.o

bheap.o : setup
	$(CC) $(CFLAGS) src/profiles/bheap.c -Isrc/include -o build/bheap.o

idr.o : setup
	$(CC) $(CFLAGS) src/profiles/idr.c -Isrc/include -o build/idr.o

alignio.o : setup
	$(CC) $(CFLAGS) src/profiles/alignio.c -Isrc/include -o build/alignio.o


# Compile shared objects for testing

tcparamprof.o : setup
	$(CC) $(CFLAGS) test/src/utest/tcparamprof.c -Isrc/include -Itest/src/utest/include -I/home/apages/tools/lcut-0.3.0/include/ -o build/tcparamprof.o

# Prepare build environment

setup:
	mkdir -p build
	mkdir -p bin

.PHONY : clean
clean :
	rm -rf build
	rm -rf bin
	rm -rf test/report/*
