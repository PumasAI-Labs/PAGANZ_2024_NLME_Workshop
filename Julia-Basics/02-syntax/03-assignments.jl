# Compound expressions
## begin..end syntax

### Multiple lines (probably the best one)
x = begin
  a = 3
  b = 2
  a + b
end

x

### Single line
x = begin a = 3; b = 2; a + b end # Not as readable, but more concise

## (;) syntax
x = (a = 3; b = 2; a + b)

# Scopes
## Global scope (accesible everywhere)
c = 2

for i = 1:5
  println("Printing $c for the $(i)th time") # We can use x inside the for loop
end

j = 1
while j <= 10
  println(j)
  global j += 1 # global lets us make sure that we are modifying the loop counter
end

## Local scope
for i = 1:5
  local d = 3
  println("$c + $d = $(c + d)")
end

d # doesn't exist

# Constants (unchaging variables)
const my_pi = 3.14
my_pi = 3 # Cannot redefine them (error for different types)
my_pi = 3.1416 # Can redefine if is the same type, but you get a warning

# Style conventions
## Variables are all lowercase
var = "my variable"
Var = "my ugly variable"

## Spaces are separated with underscores (_)
myvar = 1
my_var = 1
my_very_long_named_var = 1

## Types and Modules (next lesson) start with a capital letter and use CamelCase for spaces
Int # Type for integer
int # Does not exist
String # Type for strings
string # A function to concatenate strings
AlgebraOfGraphics # A plotting module/package/library (notice the CamelCase)
