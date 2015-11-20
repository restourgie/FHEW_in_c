q = 16
Q = 64
n = 3
N = 32
u = 9
u_2 = ((u<<2)%Q)
u_4 = ((u<<4)%Q)
u_inv = 57
t = 4

#LWE secret key
S = [-1,1,1]

#first ciphertext
a0 = [9,8,12]
m0 = 0
b0 = m0 * (q/t)

for i in range(0, n):
    b0 = (b0 + S[i]*a0[i]) % q

#second ciphertext
a1 = [2,11,13]
m1 = 1
b1 = m1 * (q/t)

for i in range(0, n):
    b1 = (b1 + S[i]*a1[i]) % q

#NAND operation
a3 = [0]*n

for i in range(0, n):
    a3[i] = (0 -(a0[i]+a1[i]))%q

b3 = ((5*q/8)-(b0+b1))%q

#DEFINE PARAMETERS
Z_Q = Integers(Q)
R.<x>  = PolynomialRing(Z_Q)
I = R.ideal([x^32 + 1])
QR = R.quotient_ring(I)

#now we first encrypt our a3 in GSW
#secret key
SN = QR(3 + 2*x + x^2 + 2*x^3 + x^4 -x^5 + 2*x^6 + x^8 + 3*x^10 -x^13 + x^14 -x^16 + x^19 -2*x^21 -2*x^22 -x^23 + x^24 + x^25 + 2*x^26 + x^27 + 3*x^28 + x^29 -x^30 + 3*x^31)
#we first encrypt a3[0] with s[0]

#generate random polynomials
##GSW CIPHERTEXT 1
gsw0 = [[QR(0)]*2 for i in range(6)]

gsw0[0][0] = 63 + 43*x + 30*x^2 + 43*x^3 + 35*x^4 + 59*x^5 + 38*x^6 + 0*x^7 + 50*x^8 + 2*x^9 + 57*x^10 + 58*x^11 + 6*x^12 + 17*x^13 + 19*x^14 + 14*x^15 + 45*x^16 + 2*x^17 + 43*x^18 + 10*x^19 + 58*x^20 + 7*x^21 + 30*x^22 + 14*x^23 + 39*x^24 + 9*x^25 + 61*x^26 + 14*x^27 + 5*x^28 + 2*x^29 + 16*x^30 + 36*x^31

gsw0[1][0] = 47 + 60*x^1 + 21*x^2 + 42*x^3 + 39*x^4 + 52*x^5 + 16*x^6 + 44*x^7 + 34*x^8 + 3*x^9 + 27*x^10 + 57*x^11 + 49*x^12 + 6*x^13 + 21*x^14 + 58*x^15 + 11*x^16 + 26*x^17 + 28*x^18 + 1*x^19 + 3*x^20 + 12*x^21 + 40*x^22 + 7*x^23 + 3*x^24 + 32*x^25 + 14*x^26 + 21*x^27 + 25*x^28 + 14*x^29 + 57*x^30 + 56*x^31

gsw0[2][0] = 26 + 59*x^1 + 45*x^2 + 11*x^3 + 31*x^4 + 28*x^5 + 12*x^6 + 4*x^7 + 2*x^8 + 61*x^9 + 20*x^10 + 26*x^11 + 47*x^12 + 9*x^13 + 8*x^14 + 63*x^15 + 28*x^16 + 23*x^17 + 17*x^18 + 24*x^19 + 46*x^20 + 11*x^21 + 3*x^22 + 44*x^23 + 54*x^24 + 48*x^25 + 59*x^26 + 34*x^27 + 57*x^28 + 7*x^29 + 30*x^30 + 19*x^31

gsw0[3][0] = 32 + 29*x^1 + 47*x^2 + 16*x^3 + 55*x^4 + 9*x^5 + 42*x^6 + 39*x^7 + 35*x^8 + 34*x^9 + 58*x^10 + 33*x^11 + 7*x^12 + 0*x^13 + 50*x^14 + 45*x^15 + 20*x^16 + 4*x^17 + 57*x^18 + 24*x^19 + 4*x^20 + 49*x^21 + 36*x^22 + 44*x^23 + 49*x^24 + 58*x^25 + 4*x^26 + 0*x^27 + 52*x^28 + 45*x^29 + 3*x^30 + 52*x^31

gsw0[4][0] = 57 + 29*x^1 + 34*x^2 + 63*x^3 + 18*x^4 + 22*x^5 + 39*x^6 + 26*x^7 + 47*x^8 + 38*x^9 + 19*x^10 + 53*x^11 + 58*x^12 + 50*x^13 + 44*x^14 + 28*x^15 + 58*x^16 + 47*x^17 + 21*x^18 + 6*x^19 + 35*x^20 + 24*x^21 + 47*x^22 + 34*x^23 + 55*x^24 + 15*x^25 + 27*x^26 + 25*x^27 + 15*x^28 + 20*x^29 + 57*x^30 + 36*x^31

gsw0[5][0] = 26 + 20*x^1 + 37*x^2 + 3*x^3 + 35*x^4 + 30*x^5 + 16*x^6 + 12*x^7 + 39*x^8 + 38*x^9 + 60*x^10 + 32*x^11 + 1*x^12 + 15*x^13 + 1*x^14 + 33*x^15 + 17*x^16 + 28*x^17 + 25*x^18 + 24*x^19 + 28*x^20 + 33*x^21 + 36*x^22 + 25*x^23 + 37*x^24 + 21*x^25 + 24*x^26 + 29*x^27 + 3*x^28 + 41*x^29 + 17*x^30 + 3*x^31

for i in range(0,6):
    gsw0[i][1] = gsw0[i][0] * SN

mm = a3[0] * S[0]
mm = (((mm % q) + q) % q) * (2*N/q)
sign = 1
if mm >= N:
    mm = mm - N
    sign = -1

gsw0[0][0] = gsw0[0][0] + sign*u*x^mm
gsw0[1][1] = gsw0[1][1] + sign*u*x^mm
gsw0[2][0] = gsw0[2][0] + sign*u_2*x^mm
gsw0[3][1] = gsw0[3][1] + sign*u_2*x^mm
gsw0[4][0] = gsw0[4][0] + sign*u_4*x^mm
gsw0[5][1] = gsw0[5][1] + sign*u_4*x^mm

#GSW CIPHERTEXT 2
gsw1 = [[QR(0)]*2 for i in range(6)]

gsw1[0][0] = 33 + 43*x^1 + 58*x^2 + 18*x^3 + 32*x^4 + 0*x^5 + 0*x^6 + 32*x^7 + 41*x^8 + 47*x^9 + 0*x^10 + 28*x^11 + 63*x^12 + 30*x^13 + 11*x^14 + 7*x^15 + 33*x^16 + 59*x^17 + 3*x^18 + 42*x^19 + 21*x^20 + 33*x^21 + 25*x^22 + 16*x^23 + 36*x^24 + 23*x^25 + 28*x^26 + 0*x^27 + 32*x^28 + 46*x^29 + 46*x^30 + 19*x^31

gsw1[1][0] = 7 + 50*x^1 + 0*x^2 + 50*x^3 + 7*x^4 + 37*x^5 + 35*x^6 + 27*x^7 + 29*x^8 + 52*x^9 + 39*x^10 + 29*x^11 + 29*x^12 + 44*x^13 + 3*x^14 + 41*x^15 + 12*x^16 + 48*x^17 + 49*x^18 + 12*x^19 + 34*x^20 + 3*x^21 + 46*x^22 + 2*x^23 + 32*x^24 + 4*x^25 + 11*x^26 + 35*x^27 + 18*x^28 + 61*x^29 + 2*x^30 + 12*x^31

gsw1[2][0] = 43 + 63*x^1 + 42*x^2 + 8*x^3 + 60*x^4 + 42*x^5 + 43*x^6 + 25*x^7 + 59*x^8 + 25*x^9 + 44*x^10 + 62*x^11 + 9*x^12 + 30*x^13 + 42*x^14 + 27*x^15 + 50*x^16 + 18*x^17 + 10*x^18 + 16*x^19 + 24*x^20 + 63*x^21 + 19*x^22 + 57*x^23 + 3*x^24 + 40*x^25 + 20*x^26 + 43*x^27 + 17*x^28 + 60*x^29 + 18*x^30 + 56*x^31

gsw1[3][0] = 23 + 51*x^1 + 53*x^2 + 62*x^3 + 52*x^4 + 52*x^5 + 28*x^6 + 23*x^7 + 17*x^8 + 53*x^9 + 45*x^10 + 46*x^11 + 17*x^12 + 38*x^13 + 56*x^14 + 26*x^15 + 34*x^16 + 43*x^17 + 63*x^18 + 48*x^19 + 15*x^20 + 58*x^21 + 48*x^22 + 50*x^23 + 0*x^24 + 7*x^25 + 48*x^26 + 17*x^27 + 37*x^28 + 29*x^29 + 30*x^30 + 51*x^31

gsw1[4][0] = 13 + 60*x^1 + 11*x^2 + 9*x^3 + 23*x^4 + 26*x^5 + 26*x^6 + 54*x^7 + 41*x^8 + 51*x^9 + 43*x^10 + 28*x^11 + 45*x^12 + 44*x^13 + 19*x^14 + 30*x^15 + 54*x^16 + 62*x^17 + 37*x^18 + 1*x^19 + 27*x^20 + 55*x^21 + 26*x^22 + 57*x^23 + 21*x^24 + 37*x^25 + 50*x^26 + 41*x^27 + 19*x^28 + 51*x^29 + 24*x^30 + 32*x^31

gsw1[5][0] = 9 + 1*x^1 + 12*x^2 + 42*x^3 + 27*x^4 + 0*x^5 + 56*x^6 + 53*x^7 + 8*x^8 + 49*x^9 + 30*x^10 + 14*x^11 + 22*x^12 + 54*x^13 + 46*x^14 + 19*x^15 + 41*x^16 + 44*x^17 + 39*x^18 + 31*x^19 + 38*x^20 + 27*x^21 + 45*x^22 + 18*x^23 + 57*x^24 + 49*x^25 + 9*x^26 + 0*x^27 + 63*x^28 + 50*x^29 + 24*x^30 + 5*x^31

for i in range(0,6):
    gsw1[i][1] = gsw1[i][0] * SN

mm = a3[1] * S[1]
mm = (((mm % q) + q) % q) * (2*N/q)
sign = 1
if mm >= N:
    mm = mm - N
    sign = -1

gsw1[0][0] = gsw1[0][0] + sign*u*x^mm
gsw1[1][1] = gsw1[1][1] + sign*u*x^mm
gsw1[2][0] = gsw1[2][0] + sign*u_2*x^mm
gsw1[3][1] = gsw1[3][1] + sign*u_2*x^mm
gsw1[4][0] = gsw1[4][0] + sign*u_4*x^mm
gsw1[5][1] = gsw1[5][1] + sign*u_4*x^mm

#GSW CIPHERTEXT 3
gsw2 = [[QR(0)]*2 for i in range(6)]

gsw2[0][0] = 42 + 33*x^1 + 5*x^2 + 11*x^3 + 4*x^4 + 59*x^5 + 43*x^6 + 4*x^7 + 42*x^8 + 1*x^9 + 59*x^10 + 23*x^11 + 15*x^12 + 24*x^13 + 43*x^14 + 5*x^15 + 28*x^16 + 62*x^17 + 10*x^18 + 17*x^19 + 6*x^20 + 31*x^21 + 37*x^22 + 32*x^23 + 13*x^24 + 14*x^25 + 60*x^26 + 8*x^27 + 44*x^28 + 52*x^29 + 44*x^30 + 16*x^31

gsw2[1][0] = 55 + 16*x^1 + 11*x^2 + 49*x^3 + 1*x^4 + 43*x^5 + 62*x^6 + 50*x^7 + 20*x^8 + 53*x^9 + 0*x^10 + 35*x^11 + 22*x^12 + 20*x^13 + 57*x^14 + 40*x^15 + 14*x^16 + 12*x^17 + 8*x^18 + 16*x^19 + 9*x^20 + 57*x^21 + 42*x^22 + 50*x^23 + 3*x^24 + 52*x^25 + 19*x^26 + 63*x^27 + 49*x^28 + 39*x^29 + 45*x^30 + 45*x^31

gsw2[2][0] = 15 + 56*x^1 + 52*x^2 + 12*x^3 + 39*x^4 + 40*x^5 + 28*x^6 + 10*x^7 + 15*x^8 + 34*x^9 + 19*x^10 + 43*x^11 + 51*x^12 + 25*x^13 + 60*x^14 + 11*x^15 + 13*x^16 + 35*x^17 + 42*x^18 + 28*x^19 + 29*x^20 + 21*x^21 + 55*x^22 + 33*x^23 + 22*x^24 + 30*x^25 + 22*x^26 + 36*x^27 + 19*x^28 + 0*x^29 + 5*x^30 + 47*x^31

gsw2[3][0] = 35 + 5*x^1 + 9*x^2 + 35*x^3 + 0*x^4 + 36*x^5 + 36*x^6 + 9*x^7 + 44*x^8 + 25*x^9 + 15*x^10 + 50*x^11 + 12*x^12 + 26*x^13 + 14*x^14 + 5*x^15 + 20*x^16 + 13*x^17 + 11*x^18 + 38*x^19 + 55*x^20 + 15*x^21 + 56*x^22 + 20*x^23 + 8*x^24 + 40*x^25 + 18*x^26 + 45*x^27 + 60*x^28 + 57*x^29 + 62*x^30 + 32*x^31

gsw2[4][0] = 55 +  29*x^1 +  53*x^2 +  62*x^3 +  12*x^4 +  0*x^5 +  31*x^6 +  44*x^7 +  26*x^8 +  62*x^9 +  23*x^10 +  13*x^11 +  55*x^12 +  34*x^13 +  58*x^14 +  52*x^15 +  21*x^16 +  16*x^17 +  51*x^18 +  63*x^19 +  34*x^20 +  5*x^21 +  26*x^22 +  40*x^23 +  1*x^24 +  47*x^25 +  1*x^26 +  27*x^27 +  20*x^28 +  53*x^29 +  0*x^30 +  6*x^31

gsw2[5][0] = 35 + 52*x^1 + 6*x^2 + 34*x^3 + 20*x^4 + 37*x^5 + 35*x^6 + 62*x^7 + 0*x^8 + 54*x^9 + 39*x^10 + 17*x^11 + 39*x^12 + 63*x^13 + 16*x^14 + 43*x^15 + 9*x^16 + 57*x^17 + 8*x^18 + 32*x^19 + 63*x^20 + 62*x^21 + 33*x^22 + 50*x^23 + 62*x^24 + 4*x^25 + 4*x^26 + 14*x^27 + 62*x^28 + 51*x^29 + 8*x^30 + 1*x^31

for i in range(0,6):
    gsw2[i][1] = gsw2[i][0] * SN

mm = a3[2] * S[2]
mm = (((mm % q) + q) % q) * (2*N/q)
sign = 1
if mm >= N:
    mm = mm - N
    sign = -1

gsw2[0][0] = gsw2[0][0] + sign*u*x^mm
gsw2[1][1] = gsw2[1][1] + sign*u*x^mm
gsw2[2][0] = gsw2[2][0] + sign*u_2*x^mm
gsw2[3][1] = gsw2[3][1] + sign*u_2*x^mm
gsw2[4][0] = gsw2[4][0] + sign*u_4*x^mm
gsw2[5][1] = gsw2[5][1] + sign*u_4*x^mm


#prepare acc
mm = (b3+(q/4))%q
mm = (((mm % q) + q) % q) * (2*N/q)
sign = 1
if mm >= N:
    mm = mm - N
    sign = -1

#init 6x2 matrix
acc = [[QR(0)]*2 for i in range(6)]

acc[0][0] = sign*u * x^mm
acc[1][1] = sign*u * x^mm
acc[2][0] = sign*u_2 * x^mm
acc[3][1] = sign*u_2 * x^mm
acc[4][0] = sign*u_4 * x^mm
acc[5][1] = sign*u_4 * x^mm

#init 6x6 matrix
_acc = [[QR(0)]*6 for i in range(6)]


mask = 0b00000011


for i in range(0,6):
    for j in range(0,2):
        temp = len(vector(acc[i][j]))
        for k in range(0,temp):
            t = int((vector(acc[i][j])[k]*u_inv)%Q)
            for l in range(0,3):
                value = mask & t
                _acc[i][j+(2*l)] += value*x^k
                t = t>>2

#add first ciphertext
acc = [[QR(0)]*2 for i in range(6)]
for i in range(0,6):
    for j in range(0,2):
        for l in range(0,6):
            acc[i][j] += (_acc[i][l] * gsw0[l][j])

#init _acc again
_acc = [[QR(0)]*6 for i in range(6)]
for i in range(0,6):
    for j in range(0,2):
        temp = len(vector(acc[i][j]))
        for k in range(0,temp):
            t = int((vector(acc[i][j])[k]*u_inv)%Q)
            for l in range(0,3):
                value = mask & t
                _acc[i][j+(2*l)] += value*x^k
                t = t>>2

#add second ciphertext
acc = [[QR(0)]*2 for i in range(6)]
for i in range(0,6):
    for j in range(0,2):
        for l in range(0,6):
            acc[i][j] += (_acc[i][l] * gsw1[l][j])

#init _acc again
_acc = [[QR(0)]*6 for i in range(6)]
for i in range(0,6):
    for j in range(0,2):
        temp = len(vector(acc[i][j]))
        for k in range(0,temp):
            t = int((vector(acc[i][j])[k]*u_inv)%Q)
            for l in range(0,3):
                value = mask & t
                _acc[i][j+(2*l)] += value*x^k
                t = t>>2

#add third ciphertext
acc = [[QR(0)]*2 for i in range(6)]
for i in range(0,6):
    for j in range(0,2):
        for l in range(0,6):
            acc[i][j] += (_acc[i][l] * gsw2[l][j])

t_TestMSB = QR(-1 + x^1 + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + x^10 + x^11 + x^12 + x^13 + x^14 + x^15 + x^16 + x^17 + x^18 + x^19 + x^20 + x^21 + x^22 + x^23 + x^24 + x^25 + x^26 + x^27 + x^28 + x^29 + x^30 + x^31)

temp_ModQ = t_TestMSB * acc[1][0]
temp_poly = t_TestMSB * acc[1][1]
temp_b = u + vector(temp_poly)[k]

#Generate SwitchingKey
switch_a = [[0]*n for i in range(N)]
switch_b = [0]*N

for i in range(0,N):
    temp_sec = vector(SN)[i]
    temp_a = vector(temp_ModQ)[i]
    switch_b[i] = -temp_sec*temp_a
    for j in range(0,3):
        switch_a[i][j] = randint(0,63)
        switch_b[i] = (switch_b[i] + (switch_a[i][j]*S[j]))%Q

#start keyswitch procedure
a_modQ = [0]*n
for i in range(0,N):
    temp_b = (temp_b - switch_b[i])%Q
    for j in range(0,n):
        a_modQ[j] = (a_modQ[j] - switch_a[i][j])%Q

#ModSwitch
type(a_modQ[0])
temp_b
res_a = [0]*n
Q
for i in range(0,n):
    res_a[i] = ((0.5) + ((a_modQ[i] * q) /Q) + q % Q)

res_b = ((0.5) + ((temp_b * q) /Q) + q % Q)









