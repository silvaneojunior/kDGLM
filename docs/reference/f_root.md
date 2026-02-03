# f_root

Calculates the root of a function given an initial value and a function
to calculate its derivatives.

## Usage

``` r
f_root(f, df, start, tol = 1e-08, n.max = 1000)
```

## Arguments

- f:

  function: A function that receives a vector and return a vector of the
  same size.

- df:

  function: A function that receives a vector and return the derivatives
  of f with respect to its arguments (if f returns a vector, it must be
  a matrix).

- start:

  vector: The initial value to start the algorithm.

- tol:

  numeric: The tolerance for the solution.

- n.max:

  numeric: The maximum number of iterations allowed.

## Value

A list containing:

- root vector: The solution for the system f(x)=0.

- f.root vector: The function f evaluated at the root.

- iter numeric: The number of steps taken.
