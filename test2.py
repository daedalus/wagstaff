from gmpy2 import *
from algos import *
import sys

exp = sys.argv[1]
n = eval(exp)
f = FSF(exp)
n1 = list_prod(f)
#print(n,f,n1)
assert n == n1
print(f)


