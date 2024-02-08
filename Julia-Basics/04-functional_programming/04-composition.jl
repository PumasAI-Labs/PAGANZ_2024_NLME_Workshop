# Often times you want to call multiple functions, passing the results to the next one

exp(sqrt(abs(-2))) # Compute abs(-2), pass the result to sqrt and then pass the result to exp

# Function composition
my_operation = (exp ∘ sqrt ∘ abs) # ∘ => \circ<TAB>
my_operation(-2)

(exp ∘ sqrt ∘ abs)(-2) # You can also do it in one line (avoid defining the function)

## More complex example: geometric mean
x = 1:10

exp(sum(log.(x)) / length(x))

geometric_mean = (exp ∘ (i -> sum(i) / length(i)) ∘ (i -> log.(i))) # Watch out: wrap anonymous functions around parenthesis
geometric_mean(x)

using StatsBase
geomean(x) # Check our results

# Function chaining (piping)
-2 |> abs |> sqrt |> exp # Opposite order from composition

## Geometric mean
x |> (i -> log.(i)) |> (i -> sum(i) / length(i)) |> exp
x .|> log |> (i -> sum(i) / length(i)) |> exp # You can vectorize the piping operator (.|>)
