# Valgrind

## Memory leak

1. For debugging information and less optimazation, give compiler appropriate options eg.) -O1 -g
2. Run
``` bash
    valgrind --leak-check=full --show-reachable=yes --show-leak-kinds=all --track-origins=yes --verbose --suppressions="suppresions.txt" --log-file="logfile" $program $arguments
```
3. Note:
  1.  Still reachables are pseudo leaks.
  1.  If openmp is enabled and used, there might be definitely or possibly lost
      memories by `pthread_create@@GLIBC` in libpthread, by `???` in libgomp. It
      seems this is not a real problem.
  1.  See the manual at www.valgrind.org


## Profile

1. For debugging information and less optimazation, give compiler appropriate options eg.) -O1 -g
2. Run
``` bash
    valgrind --tool=callgrind $program $arguments
    qcachegrind
```
