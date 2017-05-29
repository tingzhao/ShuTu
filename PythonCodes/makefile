HEADERS = -Imylib
LIBS =  -Lmylib -lmylib -lm

all: libShuTu.so ShuTuAutoTraceOneStack ShuTuAutoTrace

libShuTu.so: libShuTu.o
	gcc -shared -o libShuTu.so libShuTu.o $(LIBS)

libShuTu.o: libShuTu.c libShuTu.h
	gcc -fPIC -c libShuTu.c $(HEADERS)

ShuTuAutoTraceOneStack: ShuTuAutoTraceOneStack.c libShuTu.o
	gcc -fPIC -c ShuTuAutoTraceOneStack.c $(HEADERS)
	gcc -o ShuTuAutoTraceOneStack ShuTuAutoTraceOneStack.o libShuTu.o $(LIBS)

ShuTuAutoTrace: ShuTuAutoTrace.c libShuTu.o
	mpicc -fPIC -c ShuTuAutoTrace.c $(HEADERS)
	mpicc -o ShuTuAutoTrace ShuTuAutoTrace.o libShuTu.o $(LIBS)

test: test.c libShuTu.o
	gcc -fPIC -c test.c $(HEADERS)
	gcc -o test test.o libNeuTu.o $(LIBS)

clean:
	rm -f *.o
	rm -f *.so
