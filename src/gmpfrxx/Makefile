
include ../../config.inc

headers = $(shell ls *.h)

all : libgmpfrxx.a

libgmpfrxx.a : gmpfrxx.o mpfr_mul_d.o
	ar rs libgmpfrxx.a gmpfrxx.o mpfr_mul_d.o

%.o : %.cc $(headers)
	gcc -c $(CFLAGS) $(INC) $<

%.o : %.cpp $(headers)
	g++ -c $(CPPFLAGS) $(INC) $<

%.o : %.c $(headers)
	gcc -c $(CFLAGS) $(INC) $<

.PHONY : clean
clean:
	rm -f *.o libgmpfrxx.a
