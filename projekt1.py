# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 16:02:09 2023

@author: dawid
"""
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

    print(wyniki)
    
   