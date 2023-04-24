import numpy as np
from math import *

def hirvonen(X,Y,Z,a,e2):
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
    return(f,l,h)

def flh2XYZ(f,l,h,a,e2):
    Np = a / np.sqrt(1 - e2*np.sin(f)**2)
    N = Np (f,a,e2)
    X = (N + h) * cos(f) * cos(l)
    Y = (N + h) * cos(f) * sin(l)
    Z = (N * (1 - e2) + h) * sin(f)
    return(X,Y,Z)

def Rneu(f,l):
    R = np.array([[-np.sin(f) * np.cos(l), -np.sin(l), np.cos(f) * np.cos(l)],
                  [-np.sin(f) * np.sin(l),  np.cos(l), np.cos(f) * np.sin(l)],
                  [np.cos(f), 0., np.sin(f)]])
    return(R)
def XYZ2neu(dX,f,l):
    R = Rneu(f,l)
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

def sig(f,a,e2):
    A0 = 1 - e2/4 - 3 * e2**2/64 - 5 * e2**3/256
    A2 = 3/8 * (e2 + e2**2/4 + 15 * e2**3/128)
    A4 = 15/256 * (e2**2 + 3 * e2**3/4)
    A6 = 35 * e2**3/3072
    si = a * (A0 * f - A2 * np.sin(2 * f) + A4 * np.sin(4 * f) - A6 * np.sin(6 * f))
    return (si)
def fl2xgkygk(f,l,l0,a,e2):
    b2 = a**2 * (1-e2)
    er2 = (a**2 - b2)/b2
    dl = l - l0
    t = np.tan(f)
    eta2 = er2 * np.cos(f)**2
    Np = a / np.sqrt(1 - e2*np.sin(f)**2)
    N = Np(f,a,e2)
    xgk = sig(f,a,e2) + dl**2/2 * N * np.sin(f) * np.cos(f) * (1 + dl**2/12 * np.cos(f)**2 * (5 - t**2 + 9 * eta2 + 4 * eta2**2) + dl**4/360 * np.cos(f)**4 * (5 - 18 * t**2 + t**4 + 14 * eta2 - 330 * eta2))
    ygk = dl * N * np.cos(f) * (1 + dl**2/6 * np.cos(f)**2 * (1 - t**2 + eta2) + dl**4/120 * np.cos(f)**4 * (5 - 18 * t**2 + t**4 + 14* eta2 - 58 * eta2 * t**2))
    return(xgk, ygk)
def Tgk2000(xgk,ygk,n):
    m2000 = 0.999923
    x2000 = xgk * m2000
    y2000 = ygk * m2000 + n * 1000000 + 500000
    return(x2000,y2000)

def Tgk1992(xgk,ygk):
    m1992 = 0.9993
    x1992 = xgk * m1992 - 5300000
    y1992 = ygk * m1992 + 500000
    return(x1992,y1992)

def T2000gk(x2000,y2000,n):
    m2000 = 0.999923
    xgk = x2000/m2000
    ygk = (y2000 - n * 1000000 - 500000)/m2000
    return(xgk,ygk)
    
def T1992gk(x1992,y1992):
    m1992 = 0.9993
    xgk = (x1992 + 5300000)/m1992
    ygk = (y1992 - 500000)/m1992
    return(xgk,ygk)

#praca na projekcie