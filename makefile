CC=g++

CFLAGS=-I ./AdvXMLParser -O
LFLAGS=-L ./AdvXMLParser -ladvxml
SOURCES=experiment.cpp virtexp.cpp utils.cpp processor.cpp 
INCLUDES=virtexp.h utils.h GAEngine.h processor.h GATESTER.h 


all: experiment

experiment: $(SOURCES) $(INCLUDES)
	$(CC) $(CFLAGS) $(SOURCES) -o experiment $(LFLAGS)

clean:
	rm -f experiment *~ *.o
