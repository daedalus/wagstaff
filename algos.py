from gmpy2 import *


def find_period_1p(p):
  n = 1
  while pow(10, n, p) != 1:
    n+=1
  return n


print(find_period_1p(49))


def A002371(L):
  A=[]
  L=100
  p = 1
  while p<=L:
    p = next_prime(p)
    if p == 5 or p == 2:
      A.append(0)
    else:
      A.append(find_period_1p(p))
  return A


print(A002371(100))


def congruence(a,b,p):
  return ((a-b) % p) == 0


print(congruence(0,14,7))
print(congruence(7,14,7))
print(congruence(14,7,7))


from functools import reduce
def list_prod(lst):
  return reduce(lambda x, y: x * y, lst)


def fastpowmod(n, e, m):
  y = 1
  z = n
  while e > 0:
    if e & 1 == 1:
      y *= z % m
    z = (z*z) % m
    e >>= 1
  return y % m


print(fastpowmod(4,14,5))
print(pow(4,14,5))


def _invert(n, p):
  return fastpowmod(n, p - 2, p)


def CRT(c,m):
  M = list_prod(m)
  S = 0
  for i in range(0,len(c)):
    ni = M // m[i]
    bi = _invert(ni, m[i])
    S += (ni * bi * c[i]) 
  return S % M


print(CRT([4,24],[7,103]))
print(CRT([2,3,2],[3,5,7]))


def euler_criterion_is_prime(p):
  return fastpowmod(2, p - 1, p) == 1

print(euler_criterion_is_prime(162167))

def phi(n):
  return sum(1 for i in range(1,n) if gcd(n,i) == 1)
    

def _primitiveroot(m):
  t = phi(m)
  g = 0
  return sum(1 for g in range(2, t) if fastpowmod(g, t, m) == 1)
    #print(g,t,m)
  #  g+=1
  #return g

def primitiveroot(m):
  return phi(phi(m))

m=5
print(phi(m))
print(_primitiveroot(m))
print(primitiveroot(m))


