CC = gcc
OBJS = build/xcorr.o build/sdrnaprofiles.o
OBJS_PLOT = build/gnuplot_i.o build/plotprofiles.o

all : benchmark

benchmark : compile
	./sdrnaprofiles bnm/positive_set.bed bnm/positive_set_r1.bg bnm/positive_set_r2.bg -c 0.01 -m 18 -M 30 report

test : debug
	perl test/src/unit_test.pl test/data/unit_test.cases test/report/unit_test.dat ./sdrnaprofiles test/data/positive_set.bam test/report
	cat test/report/unit_test.dat
 
compile : xcorr.o sdrnaprofiles.o gnuplot_i.o plotprofiles.o
	$(CC) -o sdrnaprofiles $(OBJS) -Llib/ -lgsl -lgslcblas -lm -lz -lpthread -lbam -lapcluster
	$(CC) -o plotprofiles $(OBJS_PLOT) -Llib/

debug : setup
	$(CC) -c src/xcorr.c -Isrc/include/ -o build/xcorr.o -DDEBUG 
	$(CC) -c src/sdrnaprofiles.c -Isrc/include/ -Isrc/include/samtools -o build/sdrnaprofiles.o -DDEBUG
	$(CC) -o sdrnaprofiles $(OBJS) -Llib/ -lgsl -lgslcblas -lm -lz -lpthread -lbam -lapcluster

plotprofiles.o: setup
	$(CC) -c src/plotprofiles.c -Isrc/include -o build/plotprofiles.o

sdrnaprofiles.o : setup
	$(CC) -c src/sdrnaprofiles.c -Isrc/include/ -Isrc/include/samtools -o build/sdrnaprofiles.o

xcorr.o : setup
	$(CC) -c src/xcorr.c -Isrc/include/ -o build/xcorr.o

gnuplot_i.o: setup
	$(CC) -c src/gnuplot_i.c  -Isrc/include  -o build/gnuplot_i.o

setup:
	mkdir -p build

.PHONY : clean
clean :
	rm -rf build/*
	rm sdrnaprofiles
	rm plotprofiles
	rm -rf test/report/*
