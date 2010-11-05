
BASE_FLAGS=-fPIC -Wall -g -O3

CFLAGS=$(shell gsl-config --cflags) \
       $(shell python-config --cflags) \
       --std=c99 $(BASE_FLAGS)

CPPFLAGS=$(shell gsl-config --cflags) \
         $(shell python-config --cflags | sed s/-Wstrict-prototypes// ) \
         $(BASE_FLAGS)

LIBS=$(shell gsl-config --libs) $(shell python-config --libs) -lz -lgmp -lmpfr -lgmpfrxx -lnlopt
INC=$(shell python-config --includes)


