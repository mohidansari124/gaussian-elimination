# Gaussian elimination solver

This C++17 program reads a linear system from a `.lin` file and writes its solution to a `.sol` file beside the input. Existing solution files are overwritten on success.

## Build

With GCC (including MinGW on Windows):

```sh
g++ -std=c++17 -O2 -Wall -Wextra -pedantic -ffp-contract=off gaussian.cpp -o gaussian
```

On Windows this creates `gaussian.exe`. In PowerShell, invoke it as `.\gaussian.exe`.

## Run

```sh
gaussian sys1.lin
gaussian --spp sys1.lin
gaussian --double sys1.lin
gaussian --spp --double sys1.lin
```

The default is naive Gaussian elimination with single precision (`float`). `--spp` selects scaled partial pivoting; `--double` selects double precision (`double`). Flags can appear in either order, before or after the filename. Quote filenames containing spaces.

## File format

The input contains a positive integer `n`, followed by `n` rows of `n` coefficients, then `n` right-hand-side values. All values are separated by whitespace. For example, `sys1.lin` contains:

```text
3
2 1 -1
-3 -1 2
-2 1 2
8 -11 -3
```

Its solution, written to `sys1.sol`, is:

```text
2 3 -1
```

Output is one line of space-separated numbers, with enough significant digits to preserve the selected precision.

## Algorithm and limitations

Both methods perform forward elimination followed by back substitution. Naive elimination never exchanges rows. Scaled partial pivoting computes each original row's maximum absolute coefficient, then selects the remaining row maximizing `abs(a[i][k]) / scale[i]` at each step. Scales move with their rows and are not recomputed.

All arithmetic uses the selected type. The build command disables fused multiply-add contraction so multiplication and subtraction have separate rounding. On the tested compiler, `float` and `double` are 32-bit and 64-bit floating-point types.

Invalid input, zero pivots, and non-finite arithmetic produce an error on standard error and a nonzero exit status. No solution file is opened until the computation succeeds; an older solution file remains if computation fails. A zero pivot in naive elimination can occur even for a nonsingular system; try `--spp`. Singular systems do not have a unique solution. Floating-point rounding and ill-conditioning can affect accuracy, and tiny nonzero pivots are not rejected using an arbitrary tolerance.

The algorithm uses O(n^3) time and O(n^2) memory.

## Tests

After building, run `py test_gaussian.py` on Windows or `python3 test_gaussian.py` on other platforms. The tests use only Python's standard library.
