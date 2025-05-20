# TNK/HNK Source Code

To compile this code, you must have a working C++ compilation environment up to C++11. 

Issue command `g++ -O3 -I. cryptrand.c landscape.c mt19937.c -o landscape.exe -lm` from the folder that holds the contents of this repo to rebuild.

The different functionalities of this code are switched at compile-time by editing which function is run in `main()` within `landscape.c`. Further, when sweeping, the bounds of the sweep are updated within the sweep function itself (ie. `min_N`, `max_N`, and `min_K` must be changed in `run_sweep_test()` to change the bounds of a sweep). 

The seeding pattern used in `adaptive_test_mod_NK()` is deterministic and avoids state leakage from one set of $N$, $K$ parameters to the next. The function `generate_random_64bit_from_timestamp()` can be modified in `cryptrand.c` to create time-dependent seeds, but this is not necessary and would be non-deterministic. The current approach that consecutively calls the underlying `rand()` function is strong enough to avoid state leakage.