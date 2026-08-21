# Warn about cells with too little history to decompose

Warn about cells with too little history to decompose

## Usage

``` r
warn_decompose(short, thin, total, degree)
```

## Arguments

- short:

  cells with fewer steps than the trend needs

- thin:

  cells whose seasonal term rests on one year

- total:

  number of cells

- degree:

  the requested polynomial degree

## Value

invisible `NULL`
