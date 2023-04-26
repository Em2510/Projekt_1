
import csv
import argparse
import numpy as np
from math import *

def Np(f, a, e2):
    N = a / np.sqrt(1 - e2 * sin(f)**2)
    return(N)

def Mp(f,a,e2):
    M = a * (1-e2) / np.sqrt((1-e2 * np.sin(f)**2)**3)
    return(M)

def xyz2flh(X, Y, Z, a, e2):
    P = sqrt(X**2 + Y**2)
    f = np.arctan(Z/(P*(1 - e2)))
    
    while True:
        N = Np(f, a, e2)
        h = P / cos(f) - N
        fp = f
        f = np.arctan(Z/(P* (1 - e2 * N / (N + h))))
        if abs(fp - f) < (0.000001/206265):
            break
    
    l = np.arctan2(Y, X)
    return f, l, h

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='dane.txt')
    
    args = parser.parse_args()

    dane = np.loadtxt(args.input)
    wyniki = []
    with open(args.input, 'r') as plik:
        lines = plik.readlines()
        i = 0
        wyniki = []
        while i < len(lines):
            X, Y, Z, a, e2 = [float(x) for x in lines[i].strip().split()]
            wynik = xyz2flh(X, Y, Z,a, e2)
            wyniki.append(wynik)
            i += 1
            
    print('XYZ2FLH=',wynik)
    np.savetxt("flh.txt", wynik)

dane = dane[3:]
dane = ' '.join(str(x) for x in dane)


wynik = ' '.join(str(x) for x in wynik)
flhae2 = [wynik,dane]

flhae2 = ' '.join(str(x) for x in flhae2)
with open("flhae2.txt", "w") as plik:
    
    plik.write(flhae2)
    


def flh2xyz(f, l, h, a, e2):
    N = Np(f, a, e2)
    X = (N + h) * cos(f) * cos(l)
    Y = (N + h) * cos(f) * sin(l)
    Z = (N + h - N * e2) * sin(f)
    return(X, Y, Z)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default=flhae2)

   
    args = parser.parse_args()

    dane = np.loadtxt(args.input)
 
    with open(args.input, 'r') as plik:
        lines = plik.readlines()
        i = 0
        wyniki = []
        while i < len(lines):
            f, l, h, a, e2 = [float(x) for x in lines[i].strip().split()]
            wynik1 = flh2xyz(f, l, h, a, e2)
            wyniki.append(wynik1)
            i += 1
        print('FLH2XYZ=',wynik1)

def sigma(f, a, e2):
    A0 = 1 - e2/4 - 3 * e2**2/64 - 5 * e2**3/256
    A2 = (3/8) * (e2 + e2**2/4 + 15*e2**3/128)
    A4 = (15/256) * (e2**2 + (3 * e2**3)/4)
    A6 = 35 * e2**3/3072
    sigma = a * (A0*f - A2*sin(2*f) + A4*sin(4*f) - A6*sin(6*f))
    return sigma
la = l*180/pi

ns = 0 
if 16.5 > la > 14:
   ns = ns + 5
if 19.5 > la > 16.5:
   ns = ns + 6
if 22.5 > la > 19.5:
   ns = ns + 7
if 24 > la > 22.5:
   ns = ns + 8

l0 = ns/3

def fl2pl1992(f,l,a,e2,l0, m0 = 0.9993):
    b2 = a**2*(1 - e2)
    ep2 = (a**2 - b2)/b2
    dl = l - l0
    t = tan(f)
    n2 = ep2 * cos(f)**2
    N = Np(f,a,e2)
    sigm = sigma(f,a,e2)
    xgk = sigm + (dl**2/2) * N * sin(f)*cos(f)*(1 + (dl**2/12)*cos(f)**2*(5-t**2+9*n2+4*n2**2)+ ((dl**4)/360)*cos(f)**4*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*t**2))
    ygk = dl*N*cos(f)*(1+(dl**2/6)*cos(f)**2*(1 - t**2 + n2) + (dl**4/120)*cos(f)**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))
    x92 = xgk * m0 - 5300000
    y92 = ygk * m0 + 500000
    return xgk,ygk,x92,y92
xgk, ygk, x92, y92 = fl2pl1992(f, l, a, e2, l0, m0 = 0.9993)
   
print("Współrzędne X i Y w układzie 1992:", x92, y92)

def fl2pl2000(f,l,a,e2,ns ,m0= 0.999923):
    if ns == 5:
        l0 = radians(15)
    elif ns == 6:
        l0 = radians(18)
    elif ns == 7:
        l0 = radians(21)
    elif ns == 8:
        l0 = radians(24)
    b2 = a**2*(1 - e2)
    ep2 = (a**2 - b2)/b2
    dl = l - l0
    t = tan(f)
    n2 = ep2 * cos(f)**2
    N = Np(f,a,e2)
    sigm = sigma(f,a,e2)
    xgk = sigm + (dl**2/2) * N * sin(f)*cos(f)*(1 + (dl**2/12)*cos(f)**2*(5-t**2+9*n2+4*n2**2)+ ((dl**4)/360)*cos(f)**4*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*t**2))
    ygk = dl*N*cos(f)*(1+(dl**2/6)*cos(f)**2*(1 - t**2 + n2) + (dl**4/120)*cos(f)**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))
    x2000 = xgk * m0
    y2000 = ygk * m0 + ns * 1000000 + 500000
    return xgk,ygk,x2000,y2000

xgk, ygk, x20, y20 = fl2pl2000(f,l,a,e2,ns ,m0= 0.999923)
   
print("Współrzędne X i Y w układzie 2000:", x20, y20)
        