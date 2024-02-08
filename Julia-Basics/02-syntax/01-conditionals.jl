# Conditional statements allow controlling the program's flow of execution
# based on the value of a boolean expression

x = 1
y = 2

# Created with the if-elseif-else keywords
if x < y # Checks if x is less than y
  println("x is less than y")
elseif x > y # If x is not less than y, checks if it is greater than y
  println("x is greater than y")
else # If none of the above are true
  println("x is equal to y")
end # Watch out: don't forget to add the end keyword

## Try it again
x = 3
x = 2

## Tip: you can include multiple elseif
a = 12.5

if a < 5
  message = "a is less than 5" # You can do anything inside the if-elseif-else blocks
elseif a < 10 # Equivalent to 5 < x < 10
  message = "a is less than 10, but greater than 5"
elseif a < 15
  message = "a is less than 15, but greater than 10"
  # You don't necessarily have to include else (or elseif)
end

println(message)

# Ternary operator (<boolean> ? <run if true> : <run if false>)
## Convenient way to create conditional statements in a single line (if-else only)
b = 8
c = b < 10 ? 2b : b # If b is less than 10, return 2*b, otherwise just return b

## Try it again
b = 12
