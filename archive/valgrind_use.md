# Valgrind

## Memory leak

1. To get debugging information with less optimazation, give compiler appropriate
   options eg. `-O1 -p -ggdb2`
2. Run
   ``` bash
   valgrind \
       --leak-check=full --show-reachable=yes --show-leak-kinds=all \
       --track-origins=yes --verbose \
       --suppressions="suppressions.txt" --log-file="logfile" \
       $program $arguments
   ```
* Note
  1.  Still reachables are pseudo leaks.
  1.  If openmp is enabled and used, there might be definitely or possibly lost
      memories by `pthread_create@@GLIBC` in libpthread, by `???` in libgomp. It
      seems this is not a real problem.
  1.  See `valgrind --help`, and the manual at <https://valgrind.org>

### [Generate suppressions](https://stackoverflow.com/questions/36893494)

> They're 'leaks' caused by the startup code that can't be fixed by you and
> probably won't be fixed by the developers of the C++ runtime.  This is an
> extensive problem on Mac OS X. There are typically many allocations and a few
> tens of kilobytes of memory that are allocated by the startup code.

Write a minimal `main.cpp` so that valgrind can detect "'leaks' by the startup code".
``` cpp
int main() { return 0; }
```

Compile, run valgrind with generate suppressions option
``` bash
gcc -std=c++11 -O1 -p -ggdb2 -c main.cpp -o main
valgrind \
    --gen-suppressions=all --leak-check=full --show-leak-kinds=all \
    main 2>./suppressions.txt
```

Delete lines in `suppressions.txt` starting with `==` and `--`. In vim, `:g/^==/d |
g/^--/d`. Name all suppressions `<insert_a_suppression_name_hear>` to eg.
`macos-10.15.5_gcc-7.5.0-c++_valgrind-3.16.0.GIT_suppressions_1`


## Profile

1. To get debugging information with less optimazation, give compiler appropriate
   options eg. `-O1 -p -ggdb2`
2. Run
   ``` bash
   valgrind --tool=callgrind $program $arguments
   qcachegrind
   ```
