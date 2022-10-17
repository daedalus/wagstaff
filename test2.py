from gmpy2 import *
from algos import *
import sys

n = eval(sys.argv[1])
f= factor_special_forms(n)
n1 = list_prod(f)
assert n == n1
print(f)


