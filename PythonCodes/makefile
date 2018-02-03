HEADERS = -Imylib
LIBS =  -Lmylib -lmylib -lm
RPATH = -Wl,-rpath,miniconda/lib

all: libShuTu.so ShuTuAutoTraceOneStack ShuTuAutoTrace

libShuTu.so: libShuTu.o
	gcc -shared -o libShuTu.so libShuTu.o $(LIBS)

libShuTu.o: libShuTu.c libShuTu.h
	gcc -fPIC -c libShuTu.c $(HEADERS)

ShuTuAutoTraceOneStack: ShuTuAutoTraceOneStack.c libShuTu.o
	gcc -fPIC -c ShuTuAutoTraceOneStack.c $(HEADERS)
	gcc $(RPATH) -o ShuTuAutoTraceOneStack ShuTuAutoTraceOneStack.o libShuTu.o $(LIBS)

ShuTuAutoTrace: ShuTuAutoTrace.c libShuTu.o
	miniconda/bin/mpicc -fPIC -c ShuTuAutoTrace.c $(HEADERS)
	miniconda/bin/mpicc $(RPATH) -o ShuTuAutoTrace ShuTuAutoTrace.o libShuTu.o $(LIBS)

test: test.c libShuTu.o
	gcc -fPIC -c test.c $(HEADERS)
	gcc $(RPATH) -o test test.o libNeuTu.o $(LIBS)

clean:
	rm -f *.o
	rm -f *.so
	rm -f ShuTuAutoTraceOneStack
	rm -f ShuTuAutoTrace
	rm -f test
