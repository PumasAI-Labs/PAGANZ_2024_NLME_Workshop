include("01-syntax.jl")

# Write documentation for your functions (docstrings)
"""
This function takes a vector of time points and their corresponding
observations and returns the terminal slope calculated using the
last and the third to last points
"""
function terminal_slope(times, observations)

  dy = observations[end] - observations[end-2] # We are using the last and the second to last points to calculate the slope
  dt = times[end] - times[end-2]

  return dy / dt

end

terminal_slope(times, observations) # Behaves as you would expect

## Now try ?terminal_slope in the REPL

# Argument and return types
"""
Takes a vector of time points and their corresponding
observations and returns the AUC calculated through numerical
integration using the [trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule)
"""
function AUC(times::Vector, observations::Vector)

  auc = 0

  for i = 1:length(times)-1

    auc += (observations[i] + observations[i+1]) * (times[i+1] - times[i]) / 2

  end

  return auc

end

AUC(times, observations) # Works well
AUC(1, 10) # Doesn't work
typeof(AUC(times, observations)) # We get a Float64 by default, but we could enforce that

function AUC(times::Vector, observations::Vector)::Float64

  auc = 0

  for i = 1:length(times)-1

    auc += (observations[i] + observations[i+1]) * (times[i+1] - times[i]) / 2

  end

  return auc

end

AUC(times, observations)
typeof(AUC(times, observations))

# Default values
function terminal_slope(observations, times = [0, 1, 2, 4, 8, 12, 24]) # We define standard values for the time points

  dy = observations[end] - observations[end-2]
  dt = times[end] - times[end-2]

  return dy / dt

end

terminal_slope(observations) # We don't need to specify the time points because we provided a default value
terminal_slope(observations, times / 24) # We can always overwrite the default values if we want to

# Keyword arguments
## Arguments that can be specified through names instead of positions
function terminal_slope(times, observations; npoints = 2) # Now we can specify how far back we want to go to calculate the slope

  dy = observations[end] - observations[end-npoints]
  dt = times[end] - times[end-npoints]

  return dy / dt

end

terminal_slope(times, observations, 1) # Doesn't work
terminal_slope(times, observations; npoints = 1) # Keyword arguments are separated from the rest by ;
terminal_slope(times, observations, npoints = 1) # Also works, but using ; is recommended

# Anonymous functions
"""
This function takes a function and applies it to
a vector of numbers
"""
function apply(func, vector) # Watch out: you cannot use "function" as an argument name
  new_vector = [func(element) for element in vector] # Array comprehension
  return new_vector
end

## Let's change the time units to days
days(hours) = hours / 24
apply(days, times)

## Now for minutes
minutes(hours) = hours * 60
apply(minutes, times)

## In those cases, we don't really care about the name of the function.
## This is a good use case for anonymous functions
apply(i -> i / 24, times) # Syntax is arg -> f(arg)
apply(i -> i * 60, times)
