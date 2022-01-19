import numpy as num
from numpy import *
import pandas as pd
import matplotlib as plt

# stałe
a=6378137 
e2=0.00669437999013

# wspołrzędne
F_lot = 40.5298
L_lot = -3.5751
H_lot = 648.00


f = open('C:/Users/a6r14/OneDrive/Dokumenty/python/projekt/dane_2.txt').readlines()
f2 = open('C:/Users/a6r14/OneDrive/Dokumenty/python/projekt/dane_2.txt', 'r+')


for s in f:
	f2.write(s.replace('\t', ";"))
    
for s in f:    
    f2.write(s.replace(",", "."))
    

f2.close()

data = num.genfromtxt('C:/Users/a6r14/OneDrive/Dokumenty/python/projekt/dane_2.txt', delimiter=';')
tab_fi = data[:, 0:1].astype(float).squeeze()
tab_lam = data[:, 1:2].astype(float).squeeze()
tab_h = data[:, 2:3].astype(float).squeeze()


# zmiana wsp z dane.txt fi,lambda, h na wsp x, y, z
def conv_geo_xyz(fi, lam, h, a, e2):
    
    fi = num.deg2rad(fi)
    lam = num.deg2rad(lam)
    N = (a / (1 - (e2) * pow(pow(num.sin(fi), 2), 0.5)))
    x = ((N + h) * num.cos(fi) * num.cos(lam))
    y = ((N + h) * num.cos(fi) * num.sin(lam))
    z = (((N * (1 - e2) + h) * num.sin(fi)))
    return (num.array((x, y, z)))


# przeliczanie wspolrzednych fi, lam, h na wspolrzedne n, e, u
def conv_geo_neu(fi1, lam1, h1, fi2, lam2, h2):
    
    pkt_1 = conv_geo_xyz(fi1, lam1, h1, a, e2)
    pkt_2 = conv_geo_xyz(fi2, lam2, h2, a, e2)

    M = num.array(
                [
                [-num.sin(fi1) * num.cos(lam1), -num.sin(lam1), num.cos(fi1) * num.cos(lam1)],
                [-num.sin(fi1) * num.sin(lam1), num.cos(lam1), num.cos(fi1) * num.sin(lam1)],
                [num.cos(fi1), 0, num.sin(fi1)]
                ]
                )
    
    M = M.transpose()

    M2 = num.array(
                [
                [pkt_2[0] - pkt_1[0]],
                [pkt_2[1] - pkt_1[1]],
                [pkt_2[2] - pkt_1[2]]
                ]
                ).squeeze()
    n_e_u = M @ M2
    
    return (n_e_u)

# obliczanie skośnej odległosci
def skos(n, e, u):
    
    return ((pow( pow(n, 2) + pow(e, 2) + pow(u, 2) , 0.5)))

def odl_z(n, e, u):
    
    if skos(n, e, u) !=0:
        return ( u/skos(n, e, u))
    else:
        return (0)
    
def azymut(e, n):
    
    if n == 0:
        return (0)
    else:
        Az = (num.rad2deg(num.arctan(e/n)))
        if ( e > 0 and n < 0) or (e < 0 and n < 0):
            return (Az + 180)
        elif e < 0 and n > 0:
            return (Az + 360)
        else:
            return (Az)

n, e, u = conv_geo_neu(F_lot, F_lot, H_lot, tab_fi, tab_lam, tab_h)
print("skośna odległosć:", skos(n, e, u))
b = num.vectorize(azymut)
print("azymut:", b(e, n))
c = num.vectorize(odl_z)
print("zenitalna odleglość:", c(n, e, u))

#  TUTAJ JESZCZE WYKRESY !!!!!!
# Ogarnąć czy się zgadza wszystko ze str i wziąć dłuższy lot może i jak będzie git to git
    

