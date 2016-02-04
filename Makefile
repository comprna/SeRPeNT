CC = gcc
CFLAGS = -O3 -c -Wall
OBJS = build/profiles.o build/paramprof.o build/bheap.o build/idr.o build/alignio.o build/trimming.o build/xcorr.o build/iofile.o build/paramclust.o build/cluster.o build/hierarchical.o build/itvltree.o build/dtw.o build/strmap.o build/profilemap.o build/annotation.o build/dclust.o build/annotate.o build/diffproc.o build/paramdiff.o build/diffprocio.o build/npstats.o

all : srnap


# Build and compile executables
 
srnap : profiles.o annotate.o diffproc.o srnap.o
	$(CC) -o bin/srnap $(OBJS) build/srnap.o -Llib/ -lgsl -lgslcblas -lm -lz -lpthread -lbam 

srnap.o : setup
	$(CC) $(CFLAGS) src/srnap.c -Isrc/include -o build/srnap.o


# Compile shared objects

diffproc.o : paramdiff.o diffprocio.o npstats.o
	$(CC) $(CFLAGS) src/diffproc/diffproc.c -Isrc/include -o build/diffproc.o

npstats.o : setup
	$(CC) $(CFLAGS) src/diffproc/npstats.c -Isrc/include -o build/npstats.o

diffprocio.o : setup
	$(CC) $(CFLAGS) src/diffproc/diffprocio.c -Isrc/include -o build/diffprocio.o

paramdiff.o : setup
	$(CC) $(CFLAGS) src/diffproc/paramdiff.c -Isrc/include -o build/paramdiff.o

annotate.o : paramclust.o xcorr.o iofile.o dtw.o hierarchical.o profilemap.o annotation.o dclust.o npstats.o
	$(CC) $(CFLAGS) src/annotate/annotate.c -Isrc/include -o build/annotate.o

profilemap.o : itvltree.o
	$(CC) $(CFLAGS) src/annotate/profilemap.c -Isrc/include/ -o build/profilemap.o

hierarchical.o : cluster.o
	$(CC) $(CFLAGS) src/annotate/hierarchical.c -Isrc/include -o build/hierarchical.o

annotation.o : strmap.o
	$(CC) $(CFLAGS) src/annotate/annotation.c -Isrc/include -o build/annotation.o

dclust.o : setup
	$(CC) $(CFLAGS) src/annotate/dclust.c -Isrc/include -o build/dclust.o

strmap.o : setup
	$(CC) $(CFLAGS) src/annotate/strmap.c -Isrc/include -o build/strmap.o

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

profiles.o : paramprof.o bheap.o idr.o trimming.o alignio.o
	$(CC) $(CFLAGS) src/profiles/profiles.c -Isrc/include -o build/profiles.o

paramprof.o : setup
	$(CC) $(CFLAGS) src/profiles/paramprof.c -Isrc/include -o build/paramprof.o

trimming.o : setup
	$(CC) $(CFLAGS) src/profiles/trimming.c -Isrc/include -o build/trimming.o

bheap.o : setup
	$(CC) $(CFLAGS) src/profiles/bheap.c -Isrc/include -o build/bheap.o

idr.o : setup
	$(CC) $(CFLAGS) src/profiles/idr.c -Isrc/include -o build/idr.o

alignio.o : setup
	$(CC) $(CFLAGS) src/profiles/alignio.c -Isrc/include -o build/alignio.o


# Prepare build environment

setup:
	mkdir -p build
	mkdir -p bin

.PHONY : clean
clean :
	rm -rf build
	rm -rf bin
	rm -rf test/report/*
