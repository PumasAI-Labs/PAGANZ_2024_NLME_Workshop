# Filter elements from an array using a function

# filter returns all elements that satisfy a certain condition
x = 1:10

filter(iseven, x) # Get all even numbers (iseven(x) = true)
filter(!iseven, x) # Get all odd numbers (negate iseven)
filter(isodd, x) # Get all odd numbers (isodd(x) = true)

## Common use case: remove missing values
obs = [missing, 0.067, 110, 220, 220, missing, 110, 58, missing, 76]

filter(i -> i != missing, obs) # Doesn't work, we need to use the ismissing function
filter(!ismissing, obs) # Negate ismissing to get all the values that are not missing

# Count occurrences that satisfy a certain condition
count(iseven, x) # There are 5 even numbers between 1 and 10
count(ismissing, obs) # Count how many missing observations are there
