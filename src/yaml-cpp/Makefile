
include ../../config.inc

obj= \
	aliascontent.o \
	conversion.o \
	emitter.o \
	emitterstate.o \
	emitterutils.o \
	exp.o \
	iterator.o \
	map.o \
	node.o \
	null.o \
	ostream.o \
	parser.o \
	parserstate.o \
	regex.o \
	scalar.o \
	scanner.o \
	scanscalar.o \
	scantag.o \
	scantoken.o \
	sequence.o \
	simplekey.o \
	stream.o \
	tag.o


libyaml-cpp.a : $(obj)
	ar rs $@ $^

%.o : %.cpp
	g++ -c $(CPPFLAGS) $(INC) $<

.PHONY : clean
clean:
	rm -f *.o *.a

