# Sometimes we want to go from a container to a single value

## reduce
x = 1:10

reduce(+, x) # Sum of all numbers in x
reduce(*, x) # Product of all numbers in x

### Of course, Julia has convenience functions for those simple cases
sum(x)
prod(x)

### Custom function
times = [0, 1, 2, 4, 8, 12, 24]
obs = [0.01, 112, 224, 220, 143, 109, 57]

"""
Simple function to calculate the AUC  with reduce
using the trapezoidal rule (https://en.wikipedia.org/wiki/Trapezoidal_rule)
"""
function AUC(accumulated, i)

  trapz_area = (obs[i] + obs[i+1]) * (times[i+1] - times[i]) / 2

  return accumulated + trapz_area

end

reduce(AUC, 1:length(obs)-1; init = 0) # Need to set the accumulated value to 0 at the start

## Common use case: reduce + map
using Statistics
reduce(+, map(i -> (i - mean(x))^2, x)) / (length(x) - 1) # Sample variance of x

### Tip: You can calculate the sample variance var from Statistics
var(x)

## mapreduce -> map + reduce, but better
mapreduce(i -> (i - mean(x))^2, +, x) / (length(x) - 1) # mapreduce(<map>, <reduce>, <collection>)
