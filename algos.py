from gmpy2 import *


def find_period_1p(p):
  n = 1
  while pow(10, n, p) != 1:
    n+=1
  return n


def A002371(L):
  A=[]
  p = 1
  for _ in range(0,L):
    p = next_prime(p)
    if p in [5, 2]:
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
  return n - 1 if is_prime(n) else sum(1 for i in range(1, n) if gcd(i, n) == 1)


def primitiveroot(m):
  t = phi(m)
  return sum(1 for _ in range(2, t) if fastpowmod(0, t, m) == 1)
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
    j += (1 << i) if w == 1 else i << 1
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
  return x + t * p


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

  for _ in range(i1, j + 1, w):
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
  while i**2 <= n:
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
  return -1 if is_power(n+1) else 0

import math
def basepow(n):
  base = 2
  while base**2 <= n:
    power = log(n) / log(base)
    ipower = int(power)
    if (abs(power - ipower)) < 0.00000001 and pow(base, ipower, n) == 0:
      return base, ipower
    base += 1

def cunningham_decompose(n):
  pm1 = _cunningham_pm1(n)
  if pm1 != 0:
    base, power = basepow(n - pm1)
    return (base, power, pm1) if base ** power + pm1 == n else None

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


def factor_perfect(n):
  q = n
  i = 0
  while q % 2 == 0:
    q >>= 1
    i += 1
  p = (1 << i)
  return (p, q, i) if ((i << 1) + 1) == q else (p, q, -1)


def trial_factor(n):
  factors=[]
  tmp = n
  p = 2
  while p <= isqrt(tmp) + 1: 
    if tmp % p == 0:
      while tmp % p == 0:
        tmp //= p
        factors += [p]
      if is_prime(tmp):
        factors += [tmp]
    p = next_prime(p)
  return sorted(factors)


def factor_special_forms(n):
  """
  it will recursively factor numbers of the form b^m+-1 or x^2y-y^2x
  """
  if n == 1:
    return [1]
  if is_prime(n):
    print("Prime found:",n)
    return [n]
  if is_square(n):
    print("Square found:",n)
    i2 = isqrt(n)
    if is_prime(i2):
      return [i2,i2]
    f = factor_special_forms(i2)
    return f + f
  print("factor_special_forms",n)

  #pqi = factor_perfect(n)
  #if pqi != [] and pqi[0] > 1:
  #  print(pqi)
  #  if pqi[2] > 1:
  #    print("Form (2^(p-1)*(2^p-1))")
  #  #else:
  #  #  print("Form (2^n)*q")
  #    #return [pqi[0]] + factor_special_forms(pqi[1])
  #    return [pqi[0], pqi[1]]

  # Form x^2y-y^2x 
  pq = factor_leyland(n, L = 100)
  if pq != []:
    print("Form x^2y-y^2x")
    return factor_special_forms(pq[0]) + factor_special_forms(pq[1])
  bpm1 = cunningham_decompose(n)
  if bpm1 != None:
    base, power, pm1 = bpm1
    if pm1 == -1:
      print(base, "^", power,pm1, (power-3) % 4)
      # Form 2^p-1 and p = 3 (mod 4) ==> n//(2p+1),(2p+1) 
      #if base == 2 and is_prime(power) and (power - 3) % 4 == 0:
      if is_prime(power) and (power - 3) % 4 == 0:
        print("Form base^p-1 and p = 3 (mod 4) ==> n//(2p+1),(2p+1)")
        p2 = (power << 1) + 1
        if is_prime(p2) and n % p2 == 0:
          return [p2] + factor_special_forms(n // p2)
        p = base - 1
        q = n // p
        return [p] + [q] if p > 1 else [q]
      else:
        if is_prime(power):
          print("form 2^p-1 where p is prime")
          #p,q = base - 1, base ** 4 + base ** 3 + base ** 2 + base + 1
          p = base - 1
          q = sum((base ** i) for i in range(0, power))
          #print(p,q)
          return factor_special_forms(p) + factor_special_forms(q) if p > 1 else [q]
        else:
          # Form b^2m − 1 = (b^m − 1)(b^m + 1)
          if power & 1 == 0:
            print("Form b^2m − 1 = (b^m − 1)(b^m + 1)")
            if power == 6:
              print("(x-1)(x+1)(x^2-x+1)(x^2+x+1)")
              b2 = base ** 2
              p,q,r,s = base - 1, base + 1, b2 + base + 1, b2 - base +1 
              return factor_special_forms(p) + factor_special_forms(q) + factor_special_forms(r) + factor_special_forms(s)
            else:
              p2 = power >> 1
              bp2 = base ** p2
              return factor_special_forms(bp2 - 1) + factor_special_forms(bp2 + 1)
          elif power > 2:
            print("Form base^n-1 with n composite and odd")
            r = n
            F = trial_factor(power)
            p = base - 1
            r //= p
            Q = []
            #print(F)
            for f in set(F):
              q = sum((base ** i) for i in range(0, f))
              print("Found prime:",q)
              r //= q
              Q.append(q)
            #print(Q,p,q)
            #if p > 1:
            return factor_special_forms(p) + Q + factor_special_forms(r)
    else:
      print(base, "^", power,"+",pm1)
      if (2 ** int(log2(power)) == power):
        print("Form base^(2^n) + 1")
        if (n - 1) % 4 != 0:
          return [n]
        print("Form 4n^4+1")
        x = (n - 1) >> 2
        y, z = iroot(x, 4)
        if not z:
          return [n]
        #print(x, y, z)
        a, b, c = (y * y) << 1, (y << 1), 1
        #print(a, b, c) 
        p, q = a + b + c, a - b + c
            #print(n, p, q)
        return factor_special_forms(p) + factor_special_forms(q) if p*q == n else [n]
      # Form 2^n + 1 where n = -2 (mod 4)
      if (base == 2) and (power + 2) % 4 == 0:
        print("Form 2^n + 1 where n = -2 (mod 4)")
        j = (power + 2) >> 2
        j2 = 1 << j 
        j221 = (1 << ((2 * j) - 1)) 
        p, q = j221 - j2 + 1, j221 + j2 + 1
        tmp = []
        while p % 5 == 0:
          p //= 5
          tmp += [5]
        while q % 5 == 0:
          q //= 5
          tmp += [5]
        return tmp + factor_special_forms(p) + factor_special_forms(q)
      # pm1 == 1
      # Form x^3+1 == > (x+1) * (x^2 - x + 1)
      #if power == 3:
      #  p, q = base + 1, base ** 2 - base + 1
      #  return factor_special_forms(p) + factor_special_forms(q)
      if is_prime(power):
        print("base^p+1 where p is prime")
        p = base + 1
        q = sum(((-1) ** i) *  (base ** i) for i in range(0, power))
        return factor_special_forms(p) + factor_special_forms(q) if p > 1 else (p, q)
      else:
        if power & 1 == 1:
          print("Form base^n + 1 where n is composite and odd")
          p = base + 1
        else:
          print("Form base^n + 1 where n is composite and even")
          _p = sum(1 for f in trial_factor(power) if f == 2)
          #print("p",_p)
          p = base ** (2**_p) + 1
                  #return [p,q]

        q = n // p
        return factor_special_forms(p) + factor_special_forms(q)
  pqi = factor_perfect(n)
  if pqi == [] or pqi[0] <= 1:
    return [n]
  #print(pqi)
  if pqi[2] > 1:
    print("Form (2^(p-1)*(2^p-1))")
    return [pqi[0], pqi[1]]
  else:
    print("Form (2^n)*q")
    return [pqi[0]] + factor_special_forms(pqi[1])
    #return [pqi[0], pqi[1]] 


def parse_exp(exp):
  if exp.find("-") > -1:
    x = exp.split("-")
    try:
      A,B = eval(x[0]),eval(x[1])
    except:
      A,B = 0,0
    if is_square(A) and is_square(B):
      a,b = isqrt(A), isqrt(B)
      return int(abs(a - b)), int(abs(a + b))  


def pollard_rho(n, seed=2, p=2, c=1, limit = 1000):
  f = lambda x: powmod(x, p, n) + c
  x, y, d = seed, seed, 1
  C = 0
  while d == 1 and C <= limit:
    x = f(x)
    y = f(f(y))
    d = gcd((x - y), n)
    if n > d > 1:
      return d
    C += 1
  return n


def FSF(exp):
  c = 10
  n = int(eval(exp))
  F = parse_exp(exp)
  return (F if F != None and
          (F[0] * F[1]) == n else factor_special_forms(eval(exp)))

def postprocess_factors(F):
  for f in F:
    if is_prime(f):
      F2.append(f)
    else:
      Q.append(f)
  last_f = 0
  while len(Q) > 0:
    f = Q.pop()
    if log10(f) < c:
      print("trial facoring: %d cutoff %d" % (f,c))
      tmp = trial_factor(f)
      for ft in tmp:
        if is_prime(ft):
          print("Found prime:",ft)
          F2 += [ft]
        else:
          Q += [ft]
    else:
      print("pollard_rho: %d" % f)
      p = pollard_rho(f,seed=random.randint(2, 1 << int(log2(f))),limit=100000)
      #if last_p == p: 
      #  F2 += [f]
      #  break
      q = f // p
      #last_p = p
      print("Found factor: %d" % p)
      if q > 1:
        if is_prime(p):
          print("Found prime:", p)
          F2 += [p]
        else:
          Q += [p]
        if is_prime(q):
          print("Found prime: ",p)
          F2 += [q]
        else:
          Q += [q]
      else:
      #if last_f == f:
        F2 += [f]
      #last_f = f
  return [int(p) for p in sorted(F2)]

from gmpy2 import *
import random

def is_carmichael(n, k = 50):
  if is_prime(n):
    return False
  for _ in range(0,k):
    a = random.randint(2, n - 1)
    if gcd(a,n) == 1:
      if not is_fermat_prp(n, a):
        return False
  return True


def A002997(L):
  return [i for i in range(3, L) if is_carmichael(i)]



