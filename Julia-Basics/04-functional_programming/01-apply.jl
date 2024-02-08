# Functional programming is about creating your programs by applying and composing functions

"""
A function that calculates the terminal slope using the last and
the third to last observations
"""
function terminal_slope(observations; time = [0, 1, 2, 4, 8, 12, 24]) # We are setting default time values

  dy = observations[end] - observations[end-2]
  dx = time[end] - time[end-2]

  return dy / dx

end

# Now let's suppose we wan't to calculate the terminal slope for a group of subjects
population = Dict(
  "SUBJ-1" => [0.01, 112, 224, 220, 143, 109, 57],
  "SUBJ-2" => [0.01, 78, 168, 148, 119, 97, 48],
  "SUBJ-3" => [0.01, 54, 100, 91, 73, 56, 32],
)

## We could do an array comprehension
obs = collect(values(population)) # Retrieve the observations from the dictionary
[terminal_slope(observations) for observations in obs]

## Vectorize the function
terminal_slope(obs) # We don't get the result we expected (7 values instead of 3)
terminal_slope.(obs) # f.(x) tells Julia to evaluate the function for each element
log.(obs[1]) # Works with any Julia function

### Tip: you can use the @. macro if you want to vectorize multiple function calls
abs(terminal_slope.(obs)) # Now abs fails, we need to vectorize it
abs.(terminal_slope.(obs))
@. abs(terminal_slope(obs)) # @. vectorized all the function calls for us

## Use the map function
map(terminal_slope, obs)
map(log, obs[1])

### Tip: anonymous functions are often used with map
map(x -> log.(x), obs) # Calculate the logarithm of all observations (not just obs[1])
map(i -> abs(terminal_slope(i)) < 5 ? "Less than 5" : "Greater than 5", obs) # Anonymous function + ternary operator

# Execute a function for each element (map, but when the results are not needed)
foreach(terminal_slope, obs) # We don't get anything back
typeof(foreach(terminal_slope, obs))

foreach(println, keys(population)) # Print subject IDs (we just want to print)

get_id_number(subject_id) = println("$(subject_id) has ID number $(last(subject_id))")
foreach(get_id_number, keys(population))
