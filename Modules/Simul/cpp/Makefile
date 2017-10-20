SWIG = swig/swig-3.0.2/swig
SWIGFLAGS = -Wextra -c++ -python -builtin -w454 -includeall
CXX = g++
WALL = -Wall
CXXFLAGS = -fPIC -O4 -g -std=c++0x -ffast-math
OBJS = complex.o pyry3d_cpp.o mutation/disorder.o mutation/mutation.o restraints/restraints.o restraints/symmetry.o logger/logger.o
OBJSALL = $(patsubst %.cpp, %.o, $(wildcard pdb_parser/*.cpp)) $(patsubst %.cpp, %.o, $(wildcard map/*.cpp)) $(OBJS)

_pyry3d_cpp.so: $(OBJSALL) pyry3d_cpp_wrap.o
	$(CXX) $(WALL) $(CXXFLAGS) -shared $^ -o _pyry3d_cpp.so

pyry3d_cpp_wrap.cpp: pyry3d_cpp.i
	@#warning issue is not interesting. As complex allocates memory we need
	@#to create equality operator. However, it is not supported by swig.
	@#It is acceptable for us, as we won't try to create any objects with it
	@#inside python.
	$(SWIG) $(SWIGFLAGS) -o $@ $<

pyry3d_cpp_wrap.o: pyry3d_cpp_wrap.cpp
	$(CXX) $(CXXFLAGS) -I/usr/include/python2.7 -I/usr/local/include/python2.7 -c $< -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJSALL) pyry3d_cpp_wrap.cpp pyry3d_cpp_wrap.cpp.backup _pyry3d_cpp.so
	rm -f __init__.pyc
	rm -f *.pyc
	rm -f pyry3d_cpp.py

.PHONY: clean
