# Macros are a metaprogramming tool used by many packages in Julia

# Macros start with @
@time
@show
@macroexpand

# Using macros

## Space syntax @macro arguments
@time 3 + 2 # shows you how long it takes Julia to calculate 3 + 2

@time begin # time multiple expressions
  3 + 2
  6 * 5
  4 + 2^6 - 10
end

@doc println # pulls documentation

## Parenthesis syntax
@time(3 + 2)
@doc(println)

## Tip: see what a macro does with @macroexpand
@macroexpand @time 3 + 2
