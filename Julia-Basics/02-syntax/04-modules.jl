# Modules allow you to access collections of code written by someone else

## using syntax
using Pumas # Loads Pumas' code
using Statistics # Loads the statistics module from Julia

mean([10, 13, 10]) # We can now acess mean, which is provided by Statistics

### Tip: import multiple modules in a single line
using Pumas, Statistics # Separate modules using commas

### Tip: import only some parts of the module
using Statistics: mean # Statistics has a lot of code, but this only loads the mean function

## import syntax
### Same as using, but you have to use ModuleName.function everywhere
import LinearAlgebra

norm([1, 1]) # Doesn't work
LinearAlgebra.norm([1, 1])

using LinearAlgebra
norm([1, 1]) # Now it works!
