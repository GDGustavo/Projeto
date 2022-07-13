# *****************************************************
# Variables to control Makefile operation
 
CC = g++
CFLAGS = -fopenmp
 
# ****************************************************
# Targets needed to bring the executable up to date
projeto_doc: param.o execute.o giv0.o iterd.o iterd_r.o iterN.o iterN_r.o mel.o mel_2.o proj.o makefile parameters.txt
	$(CC) $(CFLAGS) -o projeto_doc param.o execute.o iterd.o giv0.o iterN.o mel.o mel_2.o iterd_r.o proj.o iterN_r.o

execute.o: execute.cpp parameters.txt
	$(CC) $(CFLAGS) -c execute.cpp

giv0.o: giv0.cpp giv0.h
	$(CC) $(CFLAGS) -c giv0.cpp

param.o: param.cpp param.h
	$(CC) $(CFLAGS) -c param.cpp

iterd.o: iterd.cpp iterd.h
	$(CC) $(CFLAGS) -c iterd.cpp	

iterd_r.o: iterd_r.cpp iterd_r.h
	$(CC) $(CFLAGS) -c iterd_r.cpp	

iterN.o: iterN.cpp iterN.h 
	$(CC) $(CFLAGS) -c iterN.cpp

iterN_r.o: iterN_r.cpp iterN_r.h 
	$(CC) $(CFLAGS) -c iterN_r.cpp

mel.o: mel.cpp mel.h 
	$(CC) $(CFLAGS) -c mel.cpp	

mel_2.o: mel_2.cpp mel_2.h 
	$(CC) $(CFLAGS) -c mel_2.cpp

proj.o: proj.cpp proj.h
	$(CC) $(CFLAGS) -c proj.cpp
