# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 13:05:57 2023

@author: dawid
"""

import argparse
import numpy as np
from math import *
from scipy import *
class Transformacje:
    def __init__(self, model):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            e2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "kras":
            self.a = 6378245.0
            self.b = 6356863.01877305
        else:
            raise NotImplementedError(f"{model} model nie zaimplementowany")
        self.flat = (self.a - self.b) / self.a
        self.e = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.e2 = (2 * self.flat - self.flat ** 2) # eccentricity**2
        

    def xyz2flh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        
        ----
        Parameters:
        X, Y, Z : FLOAT
            współrzędne w układzie ortokartezjańskim,

        ----

        Returns:
        fi
            [stopnie dziesiętne] - szerokość geodezyjna
        lambda
            [stopnie dziesiętne] - długość geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, default
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        P = sqrt(X**2 + Y**2)
        f = np.arctan(Z/(P*(1 - self.e2)))
        
        while True:
            N = self.a / np.sqrt(1 - self.e2 * sin(f)**2)
            h = P / cos(f) - N
            fp = f
            f = np.arctan(Z/(P* (1 - self.e2 * N / (N + h))))
            if abs(fp - f) < (0.000001/206265):
                break
        
        l = np.arctan2(Y, X)
        if output == "dec_degree":
            return degrees(f), degrees(l), h
        elif output =="dms":
            f = self.dms((f))
            l = self.dms((f))
            return f"{f[0]:02d}:{f[1]:02d}:{f[2]:.2f}", f"{l[0]:02d}:{l[1]:02d}:{l[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")

    def dms(self, x):
        """
        Konwertuje wartość kąta w radianach na stopnie minuty i sekundy
        
        ----
        Prameters:
        X : FLOAT
            Wartość kąta w radianach.

        ----

        Returns:
        output [STR]
            Watość kąta w stoniach minutach i sekundach 
        """
        znak = ' '
        if x < 0:
            znak = '-'
            x = abs(x)
        x = x * 180/pi
        d = int(x)
        m = int((x - d)*60)
        s = (x - d - m/60) * 3600
        print(znak,d,m,s)

    def flh2xyz(self, f, l, h):
        """
        Funkcja odwrotna do algorytmu Hirvonena. Zamienia współrzędne geodezyjne (fi, lam, h)
        na wpółrzędne ortokartezjańskie (X, Y, Z)  
        ----
        Parameters:
        fi: FLOAT
            [radiany] - szerokość geodezyjna
        lam: FLOAT
            [radiany] - długość geodezyjna
        h: FLOAT
            [metry] - wysokość elipsoidalna

        ----

        Returns:
        X, Y, Z: FLOAT
            [metry] - współrzędne ortokartezjańskie
        output [STR] optional, default
        """
        N = self.a / np.sqrt(1 - self.e2 * sin(f)**2)
        X = (N + h) * cos(f) * cos(l)
        Y = (N + h) * cos(f) * sin(l)
        Z = (N + h - N * self.e2) * sin(f)
        return(X, Y, Z)

    def XYZ2neu(self, X, Y, Z, X0, Y0, Z0):
        """
        Funkcja przedstawia transformację współrzędnych geocentrycznych odbiornika (X, Y, Z)
        do układu topocentrycznego aby uzyskać współrzędne topocentryczne odbiornika.
        Zadanie wykorzystuje także współrzędne geocentryczne anteny (X0, Y0, Z0).
        ----
        Parameters:
        X, Y, Z : FLOAT
            współrzędne geocentryczne odbiornika (satelitów),
        X0, Y0, Z0 : FLOAT
            współrzędne geocentryczne anteny,
        ----

        Returns:
        n, e, u: FLOAT
            współrzędne ortokartezjańskie
        output [STR] optional, default

        """
        f, l, h = self.xyz2flh(X,Y,Z)
        R = np.array([[-np.sin(f) * np.cos(l), -np.sin(l), np.cos(f) * np.cos(l)],
                  [-np.sin(f) * np.sin(l),  np.cos(l), np.cos(f) * np.sin(l)],
                  [np.cos(f), 0., np.sin(f)]])
        XYZ = np.array([X,
                        Y,
                        Z])
        XYZ0 = np.array([X0,
                         Y0,
                         Z0])
        XT = XYZ-XYZ0
        neu = np.linalg.inv(np.transpose(R)@R)@XT
        return neu
    
    def fl2pl1992(self, f, l, l0=radians(19), m0 = 0.9993):
        """
        Funkcja przedstawia transformację współrzędnych geodezyjnych (fi, lam) 
        do układu współrzędnych 1992(X1992, Y1992), który oparty jest na odwzorowaniu Gaussa-Krugera.
        Ze współrzędnych geodezyjnych zostają wyznaczone płaskie współrzędne kartezjańskie (Xgk, Ygk). 
        Następnie współrzędne kartezjańskie przeliczone są do układu współrzędnych 1992
        ----
        Parameters:
        fi: FLOAT
            [radiany] - szerokość geodezyjna
        lam: FLOAT
            [radiany] - długość geodezyjna
        l0: FLOAT
            [radiany] - południk środkowy
        m0: FLOAT
            skala na południku środkowym 
        ----

        Returns:
        X1992, y1992: FLOAT
            Współrzędne w układzie współrzędnych 1992
        output [STR] optional, default

        +https://ewmapa.pl/dane/wytyczne_g-1.10.pdf
        """
        b2 = self.a**2*(1 - self.e2)
        ep2 = (self.a**2 - b2)/b2
        dl = l - l0
        t = tan(f)
        n2 = ep2 * cos(f)**2
        N = self.a / np.sqrt(1 - self.e2 * sin(f)**2)
        A0 = 1 - self.e2/4 - 3 * self.e2**2/64 - 5 * self.e2**3/256
        A2 = (3/8) * (self.e2 + self.e2**2/4 + 15*self.e2**3/128)
        A4 = (15/256) * (self.e2**2 + (3 * self.e2**3)/4)
        A6 = 35 * self.e2**3/3072
        sigm = self.a * (A0*f - A2*sin(2*f) + A4*sin(4*f) - A6*sin(6*f))
        xgk = sigm + (dl**2/2) * N * sin(f)*cos(f)*(1 + (dl**2/12)*cos(f)**2*(5-t**2+9*n2+4*n2**2)+ ((dl**4)/360)*cos(f)**4*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*t**2))
        ygk = dl*N*cos(f)*(1+(dl**2/6)*cos(f)**2*(1 - t**2 + n2) + (dl**4/120)*cos(f)**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))
        x92 = xgk * m0 - 5300000
        y92 = ygk * m0 + 500000
        return x92,y92

    def fl2pl2000(self, f, l, m0= 0.999923):
        """
        Funkcja przedstawia transformację współrzędnych geodezyjnych (fi, lam) 
        do układu współrzędnych 2000(X2000, Y2000), który oparty jest na odwzorowaniu Gaussa-Krugera.
        Ze współrzędnych geodezyjnych zostają wyznaczone płaskie współrzędne kartezjańskie (Xgk, Ygk). 
        Następnie współrzędne kartezjańskie przeliczone są do układu współrzędnych 2000.
        ----
        Parameters:
        fi: FLOAT
            [radiany] - szerokość geodezyjna
        lam: FLOAT
            [radiany] - długość geodezyjna
        m0: FLOAT
            skala na południku środkowym 
        ----

        Returns:
        X2000, y2000: FLOAT
            Współrzędne w układzie współrzędnych 2000
        output [STR] optional, default

        +https://gis-support.pl/baza-wiedzy-2/podstawy-gis/uklady-wspolrzednych-w-praktyce/
        """
       
        if l<radians(14.1400) and l>radians(24.1600):
            raise ValueError('Wartość l jest niepoprawna.')
        elif l>=radians(14.1400) and l<=radians(16.5000):
            l0 = radians(15)
        elif l>radians(16.5000) and l<radians(19.5000):
            l0 = radians(18)
        elif l>=radians(19.5000) and l<radians(22.5000):
            l0 = radians(21)
        elif l>=radians(22.5000) and l<=radians(24.1600):
            l0 = radians(24)
        else:
            l0 = radians(24)
            
                
            
        b2 = self.a**2*(1 - self.e2)
        ep2 = (self.a**2 - b2)/b2
        dl = l - l0
        t = tan(f)
        n2 = ep2 * cos(f)**2
        N = self.a / np.sqrt(1 - self.e2 * sin(f)**2)
        A0 = 1 - self.e2/4 - 3 * self.e2**2/64 - 5 * self.e2**3/256
        A2 = (3/8) * (self.e2 + self.e2**2/4 + 15*self.e2**3/128)
        A4 = (15/256) * (self.e2**2 + (3 * self.e2**3)/4)
        A6 = 35 * self.e2**3/3072
        sigm = self.a * (A0*f - A2*sin(2*f) + A4*sin(4*f) - A6*sin(6*f))
        xgk = sigm + (dl**2/2) * N * sin(f)*cos(f)*(1 + (dl**2/12)*cos(f)**2*(5-t**2+9*n2+4*n2**2)+ ((dl**4)/360)*cos(f)**4*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*t**2))
        ygk = dl*N*cos(f)*(1+(dl**2/6)*cos(f)**2*(1 - t**2 + n2) + (dl**4/120)*cos(f)**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))
        x2000 = xgk * m0
        y2000 = ygk * m0 + (l0/3) * 1000000 + 500000
        return x2000,y2000