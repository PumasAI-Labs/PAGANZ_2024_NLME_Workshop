# Loops allow the repeated evaluation of an expression

# For loops
numbers = 1:10 # Range (all numbers from 1 to 10)

## Print all numbers from 1 to 10
for number in numbers
  println(number)
end # Don't forget to add end

## Tip: you can include any expression inside the loop
### For + conditionals
for number in numbers
  if iseven(number)
    println("$number is even") # String interpolation
  else
    println("$number is odd")
  end
end # Double end: one for the loop and another one for the conditionals

# While loops
## Same example as before
counter = 1
while counter <= 10 # Will execute as long as our counter is less than or equal to 10
  println(counter)
  global counter = counter + 1 # Update global variable counter
end

## Can be very convenient when you don't necessarily have to go through every instance
### Example: find where someone is in a list
names = ["Peter", "Alice", "Juan", "Bob"]
friend = "Alice"

friend_index = 1
while names[friend_index] != friend
  global friend_index += 1 # Shorthand notation to add 1 to friend_index
end

println("$friend is in the position number $friend_index of the list")
names[friend_index] # We didn't have to go through the entire list

### Same result using for + conditional + break
for index in eachindex(names) # eachindex => get all indices associated to a vector
  if names[index] == friend
    println("$(names[index]) is the position number $index of the list")
    break # Exits the loop
  end
end # Avoid infinte loops

# Array comprehensions
## Powerful way to create arrays using for loops
x = [i for i = 1:10]
x2 = [2i for i = 1:10]
greetings = ["Hello, $name" for name in names]

## Tip: you can add conditional statements in comprehensions
even_numbers = [i for i = 1:10 if iseven(i)]
