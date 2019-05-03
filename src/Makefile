# Project: xQTLImp

CPP      = g++
CC       = gcc
OBJ      = xQTLImp.o utils.o impute.o load_info.o comp.o
LINKOBJ  =  xQTLImp.o utils.o impute.o load_info.o comp.o
BIN      =  xQTLImp

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o $(BIN) $(LIBS) -lz -lpthread 

utils.o: utils.cpp
	$(CPP) -c utils.cpp -o utils.o $(CXXFLAGS)

impute.o: impute.cpp
	$(CPP) -c impute.cpp -o impute.o $(CXXFLAGS)


xQTLImp.o: xQTLImp.cpp
	$(CPP) -c xQTLImp.cpp -o xQTLImp.o $(CXXFLAGS) 

load_info.o: load_info.cpp
	$(CPP) -c load_info.cpp -o load_info.o $(CXXFLAGS)
comp.o: comp.cpp
	$(CPP) -c comp.cpp -o comp.o $(CXXFLAGS)






