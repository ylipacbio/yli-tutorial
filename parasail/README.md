
##
Files under this directory shows yli's example of using parasail library.

my-parasail.cpp:
    - align two sequences using parasail SW local alignment
    - visualize alignment tracebacks
    - obtain alignment cigar string
    - visualize alignment cigar operations

menson.build, subprojects/parasail.wrap:
    - menson build file to create a c++ project for compiling my-parasail.cpp and create executable

makefile: by typing `make`
    - meson build will be created
    - menson build will be compiled using ninja 
    - compiled executable `my` will run and print results on screen
