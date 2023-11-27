from gmpy2 import *

def fermat(n):
    if (n-2) % 4 == 0:
        raise FactorizationError
    a = isqrt(n)
    c = 0
    while not is_square(a * a - n):
        a += 1
        c+=1
    b = isqrt(a * a - n)
    print(c, a, a * a, b, b * b)
    return a - b, a + b


def lehmer(n):
    """
    fermat based integer factorization
    """
    if (n-2) % 4 == 0:
        raise FactorizationError
    y = 1
    c= 0
    while not is_square(n + y**2):
        y += 1
        c+=1
    x = isqrt(n + y**2)
    print(c, x, x*x, y, y**2)
    return x - y, x + y


p=11
q=65537

print(fermat(p*q))
print(lehmer(p*q))
print(gcdext(p,q))
