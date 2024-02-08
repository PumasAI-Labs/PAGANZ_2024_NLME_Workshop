# Julia has several types that represent numbers,
# but the main ones are integers and floats

## Integers
one = 1 # Variable assignment
two = 2
A = 267

### Tip: define large numbers using _ as a separator
large_number = 1_234_567

## Floats (real numbers)
x = 3.14
y = 2.72

### Can create "integers"
one_float = 1.0
two_float = 2.0

### Tip: check a variable's type using typeof()
typeof(one)
typeof(one_float)

### Tip: scientific notation
1e3 # Equivalent to 1*10^3 or 1000
5.67e6 # Equivalent to 5.67*10^6 or 5_670_000

## Math operations
x + y # Addition
x - y # Subtraction
x + one # Can mix integers and floats
x * y # Multiplication
x / 2 # Division
x^2 # Exponentiation
(1 + 2)^3 * 2 + 1 # evaluated in PEMDAS order

### Tip: square roots with \sqrt<TAB>
sqrt(2)
√2

## Other numeric types
5 // 37 # Rational
1 + 2im # Complex
