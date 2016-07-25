# CODeM-Toolkit
A Generator for Uncertain Scalable Stochastic Multiobjective Test Problems

The code in this repository is a `C++` implementation of the CODeM toolkit that modifies test problems for multiobjective optimization to stochastic problems.

Some examples for stochastic problems that use different features of the toolkit can be found in `misc/examples/CODeMProblems.h`.

The application allows the user to evaluate an example problem with two sets of decision vectors: random solutions and Pareto optimal solutions for the underlying deterministic problem.

## Compile
The file `CODeM.pro` is configured for `Qt` framework, and can be compiled using `qmake`.

Additionally, `Makefile` can be used to directs `make` on how to compile and link the program when `qmake` is not available on the user's machine.
This `Makefile` is intended to be used on Linux systems or Windows systems via mingw.

### How to Use `Makefile`
1. Open a terminal and `cd` to the root of the project folder.
2. Type `make` in the terminal. This will compile the project and generate the executable.
3. To run a test, simply type `make test` after step 2.
4. Type `make clean` to remove all compiled objects and the executable.

### `Makefile` Customisation
1. Change the value of `TARGET` to use a user defined name for the executable.
2. Change the value of `CXX` to use a different compiler. By default the compiler is g++.
To compile the project using the `Makefile` type `make`
3. Change the command line arguments to `make test`. The test is configured to run the executable with the `--default` option. Specialised options can be used instead.

## Run
The application accepts command line arguments such as the problem instance, dimensions, number of solutions etc.
For a full list of options type `CODeM`, `CODeM -h` or `CODeM --help`.

The output from the program is formatted as `Matlab` syntax for convenient analysis of the results.
To write the results into a file instead of the console, use the option `-f FILENAME` or `--file FILENAME`.

## Citation
Please use the following citation when referring to this work in a sceintific publication:

S. Salomon, R.C. Purshouse, I. Giagkiozis, and P.J. Fleming. A Toolkit for Generating Scalable Stochastic Multiobjective Test Problems. In *Proceedings of the
2016 Annual Conference on Genetic and Evolutionary Computation*, Denver, CO, USA, 2016

```
@inproceedings{Salomon2016Toolkit,
    address = {Denver, CO, USA},
    author = {Salomon, S. and Purshouse, R.C. and Giagkiozis, I. and Fleming, P.J.},
    booktitle = {Proceedings of the 2016 Annual Conference on Genetic and Evolutionary Computation},
    isbn = {9781450342063},
    title = {{A Toolkit for Generating Scalable Stochastic Multiobjective Test Problems}},
    year = {2016}
}
```
