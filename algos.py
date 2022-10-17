from gmpy2 import *


def find_period_1p(p):
  n = 1
  while pow(10, n, p) != 1:
    n+=1
  return n


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


def congruence(a,b,p):
  return ((a-b) % p) == 0


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


def euler_criterion_is_prime(p):
  return fastpowmod(2, p - 1, p) == 1

    
def phi(n):
  """ Slow version of euler totient function """
  if is_prime(n):
    return n - 1
  return sum([1 for i in range(1, n) if gcd(i, n) == 1])


def primitiveroot(m):
  t = phi(m)
  g = 0
  return sum(1 for g in range(2, t) if fastpowmod(g, t, m) == 1)
    #print(g,t,m)
  #  g+=1
  #return g


def primitiveroots(m):
  return phi(phi(m))


def A010554(L):
  return [1,1] + [primitiveroots(i) for i in range(3, L + 1)]

    
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
      x = pow(r, (p + 1) >> 2, p)
      return x, p - x

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


def hensel_lift(x, r, p):
  t = invert(x << 1, p) * ((r - x ** 2) // p)
  y = x + t * p
  return y


def tonelli_p2(r, p):
  """ WIP """
  x = tonelli(r, p)
  y = hensel_lift(x[0], r, p)
  p2 = p * p
  return y, (p2 - y) % (p2)


def PHI(x,p):
  return (x**p -1) // (x-1)


def _isqrt(n):
  x = 2 ** (int(log2(n)) >> 1)
  y = (x + n//x) >> 1
  while y < x:
    x = y
    y = (x + n//x) >> 1
  return x


def fermat(n):
    assert (n-2) % 4 != 0
    a = _isqrt(n)
    t = (a << 1) + 1
    b2 = a*a - n
    while not is_square(b2):
        b2 += t
        t += 2
    a = (t - 1) >> 1
    b = _isqrt(b2)
    return a - b, a + b



def cuberoot(n):
  b, _ = iroot(n,3)
  return b

def harts(N):
  """ 
  Harts one line attack 
  taken from wagstaff the joy of factoring
  """
  #for i in range(1, cuberoot(N)):
  m = 2
  i = 1
  while not is_square(m):
    s = isqrt(N * i) + 1
    m = pow(s, 2, N)
    i += 1
  t = isqrt(m)
  g = gcd(s - t, N)
  return g, N // g
  

def find_prime(a, b):
  p = 0
  c = 1
  while not is_prime(p):
    p = a**b + c
    c += 1
  return p


def lehman(N, k, r):
  Nk4 = (k*N) << 2
  x = isqrt(Nk4) + 1
  u = x*x - Nk4
  j = (isqrt(N // k) - 1) // ((r + 1) << 2)
  if (x - k) & 1 == 0:
    i1 = 1
    u += (x << 1) + 1
    x += 1
  w = 2
  if (k-1) & 1 == 0:
    w = 4
    if (x - k + N) % 4 == 0:
      i1 += 2
      u += (x + 1) << 2
      x += 2
   
  for i in range(i1, j + 1, w):
    if is_square(u):
      y = isqrt(u)
      g = gcd(x-y, N) 
      return g, N // g
    if (k-1) & 1 == 0:
      u += (x+2) << 3 
      x += 4
    else:
      u += (x + 1) << 2
      x += 2
  
  
#n = 13
#while True:
#  print("n",n)
#  try:
#    print(lehman(299944727,12,n))
#    break
#  except:
#    pass
#  n += 1

import math
def _is_power(n):
  if (n == 1):
    return True
  i = 2
  while(i * i <= n):
    val = math.log(n) / math.log(i)
    if val == int(val):
      return True
    i += 1
  return False
 
def is_cunningham(n):
  return is_power(n-1) or is_power(n+1)

def A080262(L):
  return [i for i in range(3, L + 1) if is_cunningham(i)]



def A0(L):
 return [phi(i) for i in range(3, L + 1) if is_cunningham(i)]



#def is_cunningham(n):
#  return is_power(n-1) or is_power(n+1)

#def A(L):
# return [int(phi(i)) for i in range(3, L + 3 + 1) if is_cunningham(i)]
#print(A(200))

def _cunningham_pm1(n):
  if is_power(n-1):
    return 1
  if is_power(n+1):
    return -1
  return 0

import math
def basepow(n):
  base = 2
  while (base * base <= n):
    power = log(n) / log(base)
    ipower = int(power)
    if (abs(power - ipower)) < 0.00000001 and pow(base, ipower, n) == 0:
      return base, ipower
    base += 1

def cunningham_decompose(n):
  pm1 = _cunningham_pm1(n)
  if pm1 != 0:
    base, power = basepow(n - pm1)
    if base ** power + pm1 == n:
      return base, power, pm1
    else:
      return None

#print(cunningham_decompose((2 ** 58) + 1))


def factor_leyland(N, L=100):
  for i in range(2, L + 1):
    for j in range(2, L + 1):
      a = i ** j 
      b = j ** i
      P,Q = (a + b), (a - b)
      n = P * Q
      if n == N:
        return [P, Q]
  return []
    
def factor_special_forms(n):

  if is_prime(n):
    print("Prime found:",n)
    return [n]
  if n == 1:
    return [1] 

  print("factor_special_forms",n)

 
  # Form x^2y-y^2x 
  pq = factor_leyland(n, L = 100)
  if pq != []:
    print("Form x^2y-y^2x")
    return factor_special_forms(pq[0]) + factor_special_forms(pq[1])
  
  bpm1 = cunningham_decompose(n)
  if bpm1 != None:
    base, power, pm1 = bpm1
    if pm1 == -1:
      print(base, "^", power,pm1)
      # Form 2^p-1 and p = 3 (mod 4) ==> n//(2p+1),(2p+1) 
      if base == 2 and is_prime(power) and (power - 3) % 4 == 0:
        print("Form 2^p-1 and p = 3 (mod 4) ==> n//(2p+1),(2p+1)")
        p2 = (power << 1) + 1
        if is_prime(p2):
          return [p2] + factor_special_forms(n // p2)
        else:
          return [n]
      # Form b^2m − 1 = (b^m − 1)(b^m + 1)
      if power & 1 == 0:
        print("Form b^2m − 1 = (b^m − 1)(b^m + 1)")
        p2 = power >> 1
        bp2 = base ** p2
        #print("here",p2 - 1,bp2 + 1)
        return factor_special_forms(bp2 - 1) + factor_special_forms(bp2 + 1)
        #return [(1 << p2)-1, bp2 + 1]
      if is_prime(power):
        print("form 2^p-1 where p is prime")
        #p,q = base - 1, base ** 4 + base ** 3 + base ** 2 + base + 1
        p = base - 1
        q = 0
        for i in range(0, power):
          q += (base ** i)
          #print(base,i)
        #print(p,q)
        if p > 1:
          return factor_special_forms(p) + factor_special_forms(q)
        else:
          return p, q
    else:
      print(base, "^", power,"+",pm1)
      if (2 ** int(log2(power)) == power):
        print("Form base^(2^n) + 1")
        return [n]
      # Form 2^n + 1 where n = -2 (mod 4)
      if (base == 2) and (power + 2) % 4 == 0:
        print("Form 2^n + 1 where n = -2 (mod 4)")
        j = (power + 2) >> 2
        j2 = 1 << j 
        j221 = (1 << ((2 * j) - 1)) 
        return [j221 - j2 + 1, j221 + j2 + 1]
      # pm1 == 1
      # Form x^3+1 == > (x+1) * (x^2 - x + 1)
      #if power == 3:
      #  p, q = base + 1, base ** 2 - base + 1
      #  return factor_special_forms(p) + factor_special_forms(q)
      if is_prime(power):
        print("base^p+1 where p is prime")
        q = 0
        p = base + 1
        for i in range(0, power):
          q += ((-1) ** i) *  (base ** i)
      #  p,q = base + 1, base ** 4 - base ** 3 + base ** 2 - base + 1
        if p > 1:
          return factor_special_forms(p) + factor_special_forms(q)
        else:
          return p, q

      #if base == 3:
      #  j = 0
      #  x = 0
      #  while x != power and x <= n:
      #    j += 1 
      #    x = 3 * ((j << 1) - 1)
      #    print(n, j, x)
      #  j21 = (j << 1) - 1
      #  bj21 = base ** j21
      #  bj = base ** j
      #  p, q = bj21 - bj + 1,  bj21 + bj +1
      #  return [p, q]
  return [n]
  

from gmpy2 import *
import random

def is_carmichael(n, k = 50):
  if is_prime(n):
    return False
  for i in range(0,k):
    a = random.randint(2, n - 1)
    if gcd(a,n) == 1:
      if not is_fermat_prp(n, a):
        return False
  return True


def A002997(L):
  return [i for i in range(3, L) if is_carmichael(i)]



