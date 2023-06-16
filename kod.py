# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 13:07:09 2023

@author: dawid
"""


import argparse
import numpy as np
from math import *
from scipy import *
from definicje import *

import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='dane.txt')
    parser.add_argument('--model', type=str, default="wgs84", help='Model elipsoidy: "wgs84", "grs80", "kras"')
    parser.add_argument('--transform', type=str, help='Rodzaj transformacji: "xyz2flh", "flh2xyz", "xyz2neu", "fl2pl1992", "fl2pl2000"')
    args = parser.parse_args()

    transformacje = Transformacje(args.model)

    dane = np.loadtxt(args.input, skiprows=1, delimiter=',')
    if len(dane) == 0:
        print("Plik wejściowy jest pusty")
        sys.exit()

    if args.transform is None:
        print("Nie podano transformacji")
        sys.exit()

    wyniki = []

    with open(args.input, 'r') as plik:
        lines = plik.readlines()
        i = 1
        while i < len(lines):
            values = [float(x) for x in lines[i].strip().split(',')]

            if args.model == "wgs84":
                if args.transform == "xyz2flh":
                    X, Y, Z = values[:3]
                    wynik = transformacje.xyz2flh(X, Y, Z, output='dec_degree')
                    wyniki.append(wynik)
                elif args.transform == "flh2xyz":
                    f, l, h = values[5:8]
                    wynik = transformacje.flh2xyz(f, l, h)
                    wyniki.append(wynik)
                elif args.transform == "xyz2neu":
                    X, Y, Z, X0, Y0, Z0 = values[:6]
                    wynik = transformacje.XYZ2neu(X, Y, Z, X0, Y0, Z0)
                    wyniki.append(wynik)
                elif args.transform == "fl2pl1992":
                    f, l, h = values[5:8]
                    wynik = transformacje.fl2pl1992(f, l, l0=radians(19), m0=0.9993)
                    wyniki.append(wynik)
                elif args.transform == "fl2pl2000":
                    f, l, h = values[5:8]
                    wynik = transformacje.fl2pl2000(f, l, m0=0.999923)
                    wyniki.append(wynik)
            elif args.model == "grs80":
                if args.transform == "xyz2flh":
                    X, Y, Z = values[:3]
                    wynik = transformacje.xyz2flh(X, Y, Z, output='dec_degree')
                    wyniki.append(wynik)
                elif args.transform == "flh2xyz":
                    f, l, h = values[5:8]
                    wynik = transformacje.flh2xyz(f, l, h)
                    wyniki.append(wynik)
                elif args.transform == "xyz2neu":
                    X, Y, Z, X0, Y0, Z0 = values[:6]
                    wynik = transformacje.XYZ2neu(X, Y, Z, X0, Y0, Z0)
                    wyniki.append(wynik)
                elif args.transform == "fl2pl1992":
                    f, l, h = values[5:8]
                    wynik = transformacje.fl2pl1992(f, l, l0=radians(19), m0=0.9993)
                    wyniki.append(wynik)
                elif args.transform == "fl2pl2000":
                    f, l, h = values[5:8]
                    wynik = transformacje.fl2pl2000(f, l, m0=0.999923)
                    wyniki.append(wynik)
            elif args.model == "kras":
                if args.transform == "xyz2flh":
                    X, Y, Z = values[:3]
                    wynik = transformacje.xyz2flh(X, Y, Z, output='dec_degree')
                    wyniki.append(wynik)
                elif args.transform == "flh2xyz":
                    f, l, h = values[5:8]
                    wynik = transformacje.flh2xyz(f, l, h)
                    wyniki.append(wynik)
                elif args.transform == "xyz2neu":
                    X, Y, Z, X0, Y0, Z0 = values[:6]
                    wynik = transformacje.XYZ2neu(X, Y, Z, X0, Y0, Z0)
                    wyniki.append(wynik)
            i += 1

    if len(wyniki) == 0:
        print("Niepoprawny rodzaj transformacji lub modelu elipsoidy")
        sys.exit()

    print(args.transform.upper() + "=", wyniki)
    np.savetxt(args.transform + ".txt", wyniki)