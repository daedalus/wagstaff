from gmpy2 import *
from algos import *
import sys

p = 2
while p <= 4093:
  f = factor_special_forms(2**p-1)
  if len(f) > 1:
    print(p, f)
  p = next_prime(p)

