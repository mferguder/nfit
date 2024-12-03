include ../make.inc

TESTOBJ = test.o

test: $(TESTOBJ) ../src/include/mydll.h
	$(CXX) $(TESTOBJ) -lmydll -lfftw -lrfftw

testdfp: testdfpmin.o
	$(CXX) testdfpmin.o -lmydll -lfftw -lrfftw
%.o: %.cxx
	$(CXX) -c $(CXXFLAGS) $<
