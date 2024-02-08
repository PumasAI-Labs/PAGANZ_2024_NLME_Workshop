include("02-booleans.jl")

# Strings are created with double quotes
"Hello, world!"

'Hello, world!' # Watch out: single quotes are reserved for Chars
'A'
'B'

## Tip: create longer strings spanning more than one line with triple quotes
text = """
Lorem ipsum dolor sit amet, consectetur adipiscing elit.
Quisque mollis suscipit tincidunt. Morbi vulputate libero ex,
quis maximus nunc rutrum non.
"""

# String concatenation
greeting = "Hello"
name = "Jake"

## using string()
string(greeting, ", ", name)

## using * (same as multiplication)
greeting * ", " * name # Watch out: + doesn't work

## Interpolation
"$greeting, $name"

### Tip: you can place more complicated expressions inside of $()
"One plus two is equal to $(1 + 2)"

"One plus two is equal to $(abs(-4) - 1)" # Using functions (abs for absolute value)

# Pattern matching
contains("banana", "ana") # Check if "ana" is in "banana"
occursin("ana", "banana")

startswith("banana", "ban") # Check if "banana" starts with "ban"
endswith("banana", "ana") # Check if "banana" ends with "ana"

# Formatting
sample_text = "This is an example"

uppercase(sample_text)
lowercase(sample_text)
titlecase(sample_text)

replace(sample_text, "This is" => "That was", "an example" => "a comment")
