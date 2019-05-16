# Maxflow
Maxflow is a collection of various algorithms for finding 
the maximum flow in a network, including some parallel solvers. 


### Build guide
A C++17 compliant compiler is needed to build Maxflow. 
```bash
mkdir build && cd build
cmake .. && make
cd ..
```

### Example usage
The `maxflow` binary reads a max flow problem in the [DIMACS][Dimacs_descr] format
and solves it using the specified algorithm. The Edmonds–Karp algorithm is used in the example below.
```bash
./maxflow ek -f example.inp
``` 
For the list of possible algorithms and other options, use:
```bash
./maxflow --help
```

[Dimacs_descr]: http://lpsolve.sourceforge.net/5.5/DIMACS_maxf.htm