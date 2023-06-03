## Project EKF_GEOS3
### Author: José Manuel Ros Rodrigo.
#### Email: joros@unirioja.es
---
This repository contains a partial translation to C++ of the Initial Orbit 
Determination project written in MATLAB.
##### Compiling and running the project
To effectively compile and run this project you must keep the current folder 
structure, as well as using the GCC compiler. To compile the project run the 
following instruction from the project's root directory:
```bash
g++ src/*.cpp EKF_GEOS3.cpp
```
Once compiled, you may run it using the following instruction:
```
$ ./a.out
```
##### Running unit tests 
To effectively run the unit tests you must keep the current folder 
structure, as well as using the GCC compiler. To compile the tests file run the 
following instruction from the project's root directory:
```
$ g++ src/*.cpp test.cpp
```
Once compiled, you may run them using the following instruction:
```
$ ./a.out
```
