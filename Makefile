CC = gcc
LD = gcc
CFLAGS = -O3 -Wall -ffast-math -march=native -funroll-loops
RM = /bin/rm -f -R
OBJS = galsim.o graphics/graphics.o 
EXECUTABLE = galsim
LDFLAGS=-L/opt/X11/lib -lX11 -lm

N=3000
INPUT=input_data/ellipse_N_03000.gal
STEPS=100
DT=0.00001
GRAPHICS=0
REF_INPUT=ref_output_data/ellipse_N_03000_after100steps.gal


all:$(EXECUTABLE)

compare: compare_gal_files/compare_gal_files.c
	$(CC) -lm compare_gal_files/compare_gal_files.c -o compare

$(EXECUTABLE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(EXECUTABLE)
	
graphics/graphics.o: graphics/graphics.c graphics/graphics.h
	make -C graphics

galsim.o: galsim.c
	$(CC) $(CFLAGS) -c galsim.c

test-time: galsim
	./galsim $(N) $(INPUT) $(STEPS) $(DT) $(GRAPHICS)

test-output: test-time compare
	./compare $(N) result.gal $(REF_INPUT)

clean:
	$(RM) $(EXECUTABLE) $(OBJS) compare result.gal A3 tmpdir_for_checking A3.tar.gz && make -C graphics clean

pack: clean
	mkdir A3 && cp -R graphics galsim.c report.pdf best_timing.txt Makefile A3 && tar -czf A3.tar.gz A3

check: pack
	./check-A3.sh
