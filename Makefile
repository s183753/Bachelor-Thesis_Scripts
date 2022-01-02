CXX= g++
CPPFLAGS= 
CXXFLAGS=-Wall -std=c++11
LDFLAGS=
LDLIBS=
LINK.o=$(CXX) $(LDFLAGS)

.PHONY: clean
clean:
	-$(RM) *.o