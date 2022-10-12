from gmpy2 import *


def find_period_1p(p):
  n = 1
  while pow(10, n, p) != 1:
    n+=1
  return n


print(find_period_1p(49))


def A002371(L):
  A=[]
  p = 1
  for i in range(0,L):
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
  """
  Fast exponentiation modulo function
  """
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
  """
  Chinesse remainder
  """
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
  """ Slow version of euler totient function """
  return sum(1 for i in range(1,n) if gcd(n,i) == 1)
    

def _primitiveroot(m):
  t = phi(m)
  g = 0
  return sum(1 for g in range(2, t) if fastpowmod(g, t, m) == 1)
    #print(g,t,m)
  #  g+=1
  #return g


def primitiveroots(m):
  return phi(phi(m))


m=5
print(phi(m))
print(_primitiveroot(m))
print(primitiveroots(m))


def A010554(L):
  return [1,1] + [primitiveroots(i) for i in range(3, L + 1)]
    

print(A010554(100))


def Carmichael(n):
  """ Carmichael lambda function """
  if n == 1: return 1
  if n == 2: return 1
  l = int(log2(n))
  if (l > 2)  and (n == 2 ** l):
    print(n,l)
    return 2 ** (l - 2)
  return phi(n)


def A00(L):
  return [ Carmichael(i) for i in range(1, L + 1)]


print(A00(100))


def Legendre(a, p):
    '''
    The Legendre Sybmol.
    returns 1 if a is QR(p), or -1 if NR(p), or 0 if a divides p.
    '''
    if a % p == 0:
        return 0
    # Euler's Criterion
    return 1 if pow(a, (p - 1) >> 1, p) == 1 else -1


def tonelli(r,p):
    """ 
    Tonelli Modular Square Root:
    Not the most efficient implementation but mine 
    """
    assert legendre(r,p) == 1
    f = p - 1
    e = 0

    while f & 1 == 0:
        f >>= 1
        e += 1
   
    if e == 1:
      return pow(r, (p + 1) >> 2, p)

    for n in range(2, p):
        if legendre(n, p) == -1:
            break

    R = pow(r, f, p)
    N = pow(n, f, p)
    j = 0

    for i in range(1, e):
       RNp = (R * N) % p
       x = pow(RNp, j, p)
       y = pow(2, (e - i - 1), p)
       w = pow(x , y, p) 
       if w == 1:
         j += (1 << i)
       else:
         j += i << 1

    x = (pow(r, (f + 1) >> 1, p) * pow(N, j >> 1, p)) % p
    return x, p - x


def verify_tonelli(x, r, p):
  assert ((x[0]**2 - r ** ((p + 1) >> 1)) % p) == 0
  assert ((x[1]**2 - r ** ((p + 1) >> 1)) % p) == 0
  
  assert ((x[0]**2 - (r ** ((p - 1) >> 1) * r )) % p) == 0
  assert ((x[1]**2 - (r ** ((p - 1) >> 1) * r )) % p) == 0

  assert (((x[0] ** 2) - r) % p) == 0
  assert (((x[1] ** 2) - r) % p) == 0



x = tonelli(50, 73)
verify_tonelli(x, 50, 73)
print(x)
x = tonelli(5, 41)
verify_tonelli(x, 5, 41)
print(x)


def tonelli_pe(r, p):
  """ WIP """
  x = tonelli(r, p)
  t0 ≡ (2x[0])−1((r - x[0] ** 2) / p) 
  t1 ≡ (2x[1])−1((r - x[1] ** 2) / p) 

  y0 = x0 + t0 * p
  y1 = x1 + t1 * p

  return y0,y1


