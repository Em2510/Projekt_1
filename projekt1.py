
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
    return(f, l, h)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='dane.txt')
    parser.add_argument('--a', type=float, default=6378137.0)
    parser.add_argument('--e2', type=float, default=0.00669437999014)
    args = parser.parse_args()

    dane = np.loadtxt(args.input)
    wyniki = []
    with open(args.input, 'r') as plik:
        lines = plik.readlines()
        i = 0
        wyniki = []
        while i < len(lines):
            X, Y, Z = [float(x) for x in lines[i].strip().split()]
            wynik = xyz2flh(X, Y, Z, args.a, args.e2)
            wyniki.append(wynik)
            i += 1

    print('XYZ2FLH=',wyniki)
    np.savetxt("flh.txt", wyniki)



def flh2xyz(f, l, h, a, e2):
    N = Np(f, a, e2)
    X = (N + h) * cos(f) * cos(l)
    Y = (N + h) * cos(f) * sin(l)
    Z = (N + h - N * e2) * sin(f)
    return(X, Y, Z)
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='flh.txt')
    parser.add_argument('--a', type=float, default=6378137.0)
    parser.add_argument('--e2', type=float, default=0.00669437999014)
    args = parser.parse_args()

    dane = np.loadtxt(args.input)
 
    with open(args.input, 'r') as plik:
        lines = plik.readlines()
        i = 0
        wyniki = []
        while i < len(lines):
            f, l, h = [float(x) for x in lines[i].strip().split()]
            wynik1 = flh2xyz(f, l, h, args.a, args.e2)
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


def fl2pl1992(f,l,a,e2,l0=radians(19), m0 = 0.9993):
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='flh.txt')
    parser.add_argument('--a', type=float, default=6378137.0)
    parser.add_argument('--e2', type=float, default=0.00669437999014)
    args = parser.parse_args()
    FL = np.loadtxt(args.input, usecols=(0,1))
    
    dane = np.loadtxt(args.input)
    
    with open(args.input, 'r') as plik:
        lines = plik.readlines()
        i = 0
        wyniki = []
        while i < len(lines):
            f,l,h = [float(x) for x in lines[i].strip().split()]
            wynik = fl2pl1992(f,l,args.a, args.e2,l0=radians(19), m0 = 0.9993)
            wyniki.append(wynik[2:])
            i += 1
        print('FL2PL1992=',wyniki)   

# brakująca funkcja
def XYZ2neu(dX,X,Y,Z,a,e2):
        p = np.sqrt(X**2 + Y**2)
        f = np.arctan(Z/(p * (1 - e2)))
        Np = a / np.sqrt(1 - e2*np.sin(f)**2)
        while True:
            N = Np(f,a,e2)
            h = p / np.cos(f) - N
            fp = f
            f = np.arctan(Z / (p * (1 - e2 * N / (N + h))))
            if abs(fp - f) < (0.000001/206265):
                break
        l = np.arctan2(Y,X)
        R = np.array([[-np.sin(f) * np.cos(l), -np.sin(l), np.cos(f) * np.cos(l)],
                    [-np.sin(f) * np.sin(l),  np.cos(l), np.cos(f) * np.sin(l)],
                    [np.cos(f), 0., np.sin(f)]])
        return(R.T @ dX)
def neu2saz(dx):
        s = np.sqrt(dx @ dx)
        alfa = np.arctan2(dx[1],dx[0])
        z = np.arccos(dx[2]/s)
        return(s,alfa,z)
def saz2neu(s,alfa,z):
        dneu = np.array([s * np.sin(z) * np.cos(alfa),
                        s * np.sin(z) * np.sin(alfa),
                        s * np.cos(z)])
        return(dneu)

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
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='flh.txt')
    parser.add_argument('--a', type=float, default=6378137.0)
    parser.add_argument('--e2', type=float, default=0.00669437999014)
    args = parser.parse_args()
    FL = np.loadtxt(args.input, usecols=(0,1))
    dane = np.loadtxt(args.input)
    with open(args.input, 'r') as plik:
        lines = plik.readlines()
        i = 0
        wyniki = []
        while i < len(lines):
            f,l,h = [float(x) for x in lines[i].strip().split()]
            wynik = fl2pl2000(f,l,args.a, args.e2, ns = 5, m0 = 0.999923)
            wyniki.append(wynik[2:])
            i += 1
        print('FL2PL2000=',wyniki)    
        
        
        
