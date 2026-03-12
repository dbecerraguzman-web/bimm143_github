# Class 6: R functions
Dariana Becerra Guzman (A17506182)

- [Background](#background)
- [Our First Function](#our-first-function)
- [A second function](#a-second-function)
- [A Protein generating function](#a-protein-generating-function)

## Background

All functions in R have at least 3 things:

- A **name** that we use to call the function.
- One or more input **arguments**.
- The **body** the lines of R code that do the work.

## Our First Function

Let’s write a silly wee function called `add()` to add some numbers (the
input arguments).

``` r
add <- function(x,y) {
  x + y
}
```

Now we can use this function

``` r
add(100,1)
```

    [1] 101

``` r
add (c(100, 1, 100), 1)
```

    [1] 101   2 101

> Q. What if I give multiple element vector to `x` and `y`?

``` r
add(c(100, 1), c(100,1))
```

    [1] 200   2

> What if I give three inputs to the functions?

``` r
#add(x=c(100, 1), y=1, z=1)
```

> Q. What if I give only one input to the add function?

``` r
addnew <- function(x, y=1) {
  x + y
}
```

``` r
addnew(x=100)
```

    [1] 101

``` r
addnew(c(100,1), 100)
```

    [1] 200 101

If we write our function with input arguments having no default value
then the use will be required to set them when they use the function. We
can give out input arguments “default” values by setting the equal to
some sensible value

## A second function

Let’s try something more interesting:Make a sequence generating tool…

The `sample()` function can be a useful starting point here:

``` r
sample(1:10, size=4)
```

    [1] 4 5 8 7

> Q. Generate 9 random numbers taking from the input vector x=1:10?

``` r
sample(1:10, size=9)
```

    [1]  5  2  3  7 10  1  8  9  6

> Q. Generate 12 random numbers taking from the input vector x=1:10?

``` r
sample(1:10, size=12, replace=TRUE)
```

     [1]  4  4  4  6  1  7  8  8 10  4  2  5

> Q. Write code for the `sample()` function that generated nucleotide
> sequences of length 6?

``` r
sample(x=c("A", "G","T", "C"), size=6, replace=TRUE)
```

    [1] "A" "A" "C" "G" "T" "G"

> Q. Write a first function `generate_dna()` that returns a **user
> specified length** DNA sequence.

``` r
generate_dna <- function(len) {
  sample(x=c("A", "G", "T", "C"), size=len, replace=TRUE)
}
```

``` r
generate_dna(len=10)
```

     [1] "A" "C" "T" "A" "G" "C" "A" "G" "A" "T"

> **Key Points** Every function in R looks fundamentally the same in
> terms of its structure. Basically 3 things: name, input, and body.

    name <- function(input) {

    }

> Functions can have multiple inputs. These can be **required**
> arguments or **optional** arguments. With optional arguments having a
> set default value.

> Q. Modify and improve our `generate_dna()` function to return it’s
> generated sequence in a more standard format like “AGTATA” rather than
> the vector “A”, “C”, “G”, “T”.

``` r
generate_dna <- function(len=6, fasta=TRUE) {
  
 ans <- sample(c("A", "G", "T", "C"),
                size=len, replace=TRUE)
 if(fasta) {
   cat("Single-element vector output")
 ans <- paste(ans, collapse = "")
 } else {
   cat("Multi-element vector ouput")
 }
 
 return(ans)
}

generate_dna()
```

    Single-element vector output

    [1] "GACGCG"

``` r
generate_dna(fasta=TRUE)
```

    Single-element vector output

    [1] "ACTAGA"

The `paste()` function - it’s job is to join up or stick together
(a.k.a. paste) input strings together.

``` r
paste(c("alice", "barry"), "loves R")
```

    [1] "alice loves R" "barry loves R"

Flow control means where the R brain goes in your code

``` r
good_mood <- FALSE

if(good_mood) {
  cat("Great!")
} else {
  cat("Bummer!")
}
```

    Bummer!

## A Protein generating function

> Q. Write a function, called `generate_protein()`, that generates a use
> specified length protein sequence.

> Q. Use that function to generate randomm protein sequences between
> length 6 and 12.

> Q. Are any of your sequences unique i.e. not found anywhere in nature?

There are 20 natural amino acids

``` r
aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
```

``` r
generate_protein <- function(len, fasta=TRUE) {
 
   # The amino-acids to sample from
  aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  
   # Draw n=len amino acids to make our sequence
 ans <- sample((aa), size=len, replace=TRUE)
 ans <- paste(ans, collapse="")
 return(ans)
}
```

``` r
my_seq <- generate_protein(42)
my_seq
```

    [1] "PMKKKKMPTYSPGNNTYELERIANGHKNRMRKGTRWAADGGC"

> Q. Use that function to generate randomm protein sequences between
> length 6 and 12.

``` r
for(i in 6:12) {
  #FASTA ID line ">id"
  cat(">",i, sep="", "\n")
  # Our random protein sequence line
  cat(generate_protein(i), "\n")
}
```

    >6
    VGRQSH 
    >7
    KAKWTDW 
    >8
    NTFYRHYS 
    >9
    SDSGTVHYH 
    >10
    LHYCEQYEYI 
    >11
    KICMHIDIGWR 
    >12
    VNWECDMQIPFH 

> Q. Are any of your sequences unique i.e. not found anywhere in nature?

9, 10, 11, 12
