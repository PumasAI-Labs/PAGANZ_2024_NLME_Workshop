# Booleans are variables that can only take two values: true or false

include("01-numeric.jl")

# Comparisons
1 < 2
x < y
y ≤ 6 # \leq<TAB>, or <=
y ≥ 6 # \geq<TAB>, or >=
one == one_float # Equals (different from =, which is for assignment)
1.00000001 ≈ 1 # Approx (\approx<TAB>)

# Operators
!(1 < 2) # Negation
(1 < 2) || (2 < 1) # OR
(1 < 2) && (2 < 1) # AND

# Types
!true
true && false

# Tip: boolean values can be used for math (true = 1, false = 0)
true + true
one + true
(x > 6) * x
