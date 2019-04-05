# Project: Project_gene

CPP      = g++
CC       = gcc
OBJ      = impg.o impz.o linsubs.o main_s.o util.o zgenbt.o
LINKOBJ  = impg.o impz.o linsubs.o main_s.o util.o zgenbt.o
BIN      = Project_gene

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o $(BIN) $(LIBS) -lz -lpthread

impg.o: impg.cpp
	$(CPP) -c impg.cpp -o impg.o $(CXXFLAGS)

impz.o: impz.cpp
	$(CPP) -c impz.cpp -o impz.o $(CXXFLAGS)

linsubs.o: linsubs.cpp
	$(CPP) -c linsubs.cpp -o linsubs.o $(CXXFLAGS)

main_s.o: main_s.cpp
	$(CPP) -c main_s.cpp -o main_s.o $(CXXFLAGS)

util.o: util.cpp
	$(CPP) -c util.cpp -o util.o $(CXXFLAGS)

zgenbt.o: zgenbt.cpp
	$(CPP) -c zgenbt.cpp -o zgenbt.o $(CXXFLAGS)
