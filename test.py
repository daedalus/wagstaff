from algos import *

def test():
  print(find_period_1p(49))
  print(A002371(100))
  print(congruence(0,14,7))
  print(congruence(7,14,7))
  print(congruence(14,7,7))
  print(fastpowmod(4,14,5))
  print(pow(4,14,5))
  print(CRT([4,24],[7,103]))
  print(CRT([2,3,2],[3,5,7]))
  print(euler_criterion_is_prime(162167))
  m=5
  print(phi(m))
  print(primitiveroot(m))
  print(primitiveroots(m))
  print(A010554(100))
  print(A00(100))
  x = tonelli(50, 73)
  verify_tonelli(x, 50, 73)
  print(x)
  x = tonelli(5, 41)
  verify_tonelli(x, 5, 41)
  print(x)
  x=tonelli_p2(37987, 239)
  verify_tonelli(x, 37987, 239)
  print(x)
  print(fermat(65537*65539))
  print(harts(13290059))
  #x=     94738740796943840961823530695778701408987757287583492665919730017973847138345511139064596113422435977583856843887008168202003855906708039013487349390571801141407245039011598810542232029634564848797998534872251549660454277336502838185642937637576121533945369150901808833844341421315429263207612372324026271327
  n = find_prime(2,500) * find_prime(2,501) 
  print(n)
  print(harts(n))
  print(lehman(13290059,6,7))  
  print(A080262(100))
  print(A0(300))
  n = 2 ** 58 + 1
  f = factor_special_forms(n)
  #assert f[0]*f[1] == n
  print(n,f)

  print(is_carmichael(5))
  print(is_carmichael(15))
  print(is_carmichael(561))
  print([phi(i) for i in A002997(20000)])

#test()

def test2():
  print(factor_special_forms(2**58-1))
  print(factor_special_forms(7625597484988))
  #print(factor_special_forms(1023490369077469249537))

test2()
