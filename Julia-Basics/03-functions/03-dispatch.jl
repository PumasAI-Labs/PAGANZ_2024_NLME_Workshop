include("02-advanced.jl")

# Julia functions can behave differently depending on the types of arguments that are passed to them

# We already saw how to constrain argument types
function AUC(times::Vector, observations::Vector)

  auc = 0

  for i = 1:length(times)-1

    auc += (observations[i] + observations[i+1]) * (times[i+1] - times[i]) / 2

  end

  return auc

end

AUC(1, 10) # Throws an error because both arguments need to be Vectors
AUC(times, observations)

# Suppose we would also like to have something that calculates the AUC for a population
population = Dict(
  "SUBJ-1" => [0.01, 112, 224, 220, 143, 109, 57],
  "SUBJ-2" => [0.01, 78, 168, 148, 119, 97, 48],
  "SUBJ-3" => [0.01, 54, 100, 91, 73, 56, 32],
)

AUC(times, population) # Of course, our function won't work here

function AUC_pop(times::Vector, population::Dict)

  auc_values = Dict() # We are creating an empty dictionary to add the results

  for subject in population

    auc = 0
    observations = subject.second # Get the observations for a given subject

    for i = 1:length(times)-1

      auc += (observations[i] + observations[i+1]) * (times[i+1] - times[i]) / 2

    end

    subject_id = subject.first
    auc_values[subject_id] = auc

  end

  return auc_values

end

AUC_pop(times, population)

# But now we created another function that we need to remember
## Multiple dispatch can solve that
function AUC(times::Vector, population::Dict) # We use the same name, but different argument types

  auc_values = Dict() # We are creating an empty dictionary to add the results

  for subject in population

    auc = 0
    observations = subject.second # Get the observations for a given subject

    for i = 1:length(times)-1

      auc += (observations[i] + observations[i+1]) * (times[i+1] - times[i]) / 2

    end

    subject_id = subject.first
    auc_values[subject_id] = auc

  end

  return auc_values

end

AUC(times, population) # Now it works
AUC(times, observations) # Julia figures out which method to use depending on the arguments that we pass

AUC # Note that now "AUC" is a function with 2 methods

## Tip: check methods available for a function
methods(AUC)
methods(string) # An example of a Julia function that has many methods
