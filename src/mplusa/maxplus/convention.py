from math import prod, inf
from numpy import maximum

ADD = max
MULT = sum
POW = prod
IDENTITY = 0
ZERO = -inf
ADD_MATRICES = maximum

# Check if the value is not ZERO
def ZERO_COMPARISON(value):
    return value > ZERO

# Compare the results of operations involving addition (useful for optimization algorithms)
# Answers the question "which value is more extreme?" for the given semiring
def SOFT_EXTREME_COMPARISON(value_1, value_2):
    return value_1 >= value_2

def HARD_EXTREME_COMPARISON(value_1, value_2):
    return value_1 > value_2
