# Functions
## Powerful and commonly used tool in Julia and other programming languages

## Multiple lines syntax
function geo_mean(values) # function name(args)
  prod(values)^(1 / length(values)) # Manipulate args
end # Watch out: end is required

geo_mean(1:10)
geo_mean(rand(10)) # We have created code that works for any vector of numbers

### Tip: you can calculate geometric means with the StatsBase package
using StatsBase
geomean(1:10)

### Multiple arguments
function terminal_slope(times, observations) # Here we have two arguments, and we could add more by separating them with commas 

  dy = observations[end] - observations[end-2] # We are using the last and the second to last points to calculate the slope
  dt = times[end] - times[end-2]

  return dy / dt # The return function indicates what should be the result of calling the function

end

observations = [0.01, 112, 224, 220, 143, 109, 57]
times = [0, 1, 2, 4, 8, 12, 24]

terminal_slope(times, observations)

### No arguments
pwd() # Prints the present working directory

## Compact function assignment
geo_mean(values) = prod(values)^(1 / length(values)) # name(args) = result

geo_mean(1:10)
geo_mean(rand(10))

terminal_slope(times, observations) =
  (observations[end] - observations[end-2]) / (times[end] - times[end-2])

terminal_slope(times, observations)

## Return multiple values
using Statistics

function summary_statistics(values)

  min = minimum(values)
  max = maximum(values)
  q1 = quantile(values, 0.25)
  q2 = quantile(values, 0.5)
  q3 = quantile(values, 0.75)

  return min, max, q1, q2, q3  # Return multiple values with return1, return2, return3, ...
end

summary = summary_statistics(1:10) # We get a Tuple

# Extract values by indexing the Tuple
min = summary[1]
q2 = summary[4]

## Tip: a useful way to unpack multiple values
min, max, q1, q2, q3 = summary_statistics(1:10)
min
q2
