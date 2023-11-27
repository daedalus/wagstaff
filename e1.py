from gmpy2 import *

def express(p):
  f = 1
  e = int(log2(p)) 
  print(e) 
  while (p-1 != (2**e) * f):
    print(p-1,(2**e)*f)
    f+=2
  return f


def tonelli(r,p):
  f = p - 1
  e = 0

  while f & 1 == 0:
      f >>= 1
      e += 1

  if e == 1:
    return powmod(r, (p + 1) >> 2, p)

  for n in range(2, p):
    if legendre(n, p) == -1:
      break

  R = pow(r, f, p)
  N = pow(n, f, p)
  j = 0

  for i in range(1, e):
    x = pow((R * N), j, p)
    y = pow(2, (e - i - 1), p)
    w = pow(x , y, p)
    j += (1 << i) if w == 1 else i << 1
  x = (pow(r, (f + 1) >> 1, p) * pow(N, j >> 1, p)) % p
  return x, p - x


print(tonelli(50,73))
print(tonelli(5,41))

#print(a[0] * (2 ** a[1]) + 1)
