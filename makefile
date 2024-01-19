CFLAGS=$(shell root-config --cflags)
LIBS=$(shell root-config --libs)
OPT=-Wall

CSR_Online_Server:CSR_Online_Server.cpp
	g++ ${CFLAGS} ${LIBS} ${OPT} -o CSR_Online_Server CSR_Online_Server.cpp
