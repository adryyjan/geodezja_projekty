from random import random
import numpy as num
import random


a = 6378137
e2 = 0.0066943800290

A = [50.25, 20.75]
B = [50, 20.75]
C = [50.25, 21.25]
D = [50, 21.25] 



def srodek(x , y):
    
    return ((x + y)/2)

def dec_st(x):
    h = int(x)
    m = int((x - h) * 60)
    s = round((x - h - m/60) * 3600, 5)

    return (str(str(h) + "°"+str(m) + "'"+str(s) + "''"))

def vincent(a, e2, fi1, lam1, fi2, lam2):
    b = (a * pow((1 - e2), 0.5))
    f = (1 - b/a)
    fi1 = num.deg2rad(fi1)
    fi2 = num.deg2rad(fi2)
    lam1 = num.deg2rad(lam1)
    lam2 = num.deg2rad(lam2)
    
    L1 = lam2 - lam1
    
    U_a = num.arctan(1 - f) * num.tan(fi1)
    U_b = num.arctan(1 - f)  * num.tan(fi2)
    
    warunek = num.deg2rad(0.00001/3600)
    
    L2 = L1
    
    while abs(L1 - L2) < warunek:
        sin_sig = pow(num.cos(U_b) * pow(num.sin(L1), 2) + pow(num.cos(U_a) * num.sin(U_b) - num.sin(U_a) * num.cos(U_b) * num.cos(L1), 2), 0.5)
        
        cos_sig = num.sin(U_a) * num.sin(U_b) + num.cos(U_a) * num.cos(U_b) * num.cos(L1)
        
        sig = num.arctan(sin_sig / cos_sig)
        
        sin_alfa = num.cos(U_a) * num.cos(U_b) * num.sin(L1) / num.sin(sig)
        
        cos_al_2 = (1 - pow(sin_alfa, 2))
        
        cos_2_sig = pow(cos_al_2, 0.5 )- (2 * num.sin(U_a) * num.sin(U_b) / cos_al_2)
        
        C = f / 16 * cos_al_2 * (4 + f * (4 - 3 * cos_al_2))
        
        L2 = L1 + (1 - C) * f * sin_alfa * (sig + C * sin_sig * (cos_2_sig + C * cos_sig * (-1 + 2 * pow(cos_2_sig, 2))))
        

    L1 = L2 
    
    u_2 = (pow(a, 2) - pow(b, 2)) / pow(b, 2) * cos_al_2
    
    A = 1 + u_2 / 16384 * (4096 + u_2 * (-768 + u_2 * (320 - 175 * u_2)))
    
    B = u_2 / 1024 * (256 + u_2 * (-128 + u_2 * (74 - 47 * u_2)))
    
    del_sig = B * sin_sig * (cos_2_sig + 1 / 4 * B * (cos_sig * (-1 + 2 * pow(cos_2_sig, 2)) - 
                                                        
                                                        1 / 6 * B * cos_2_sig * (-3 + 4 * pow(sin_sig, 2)) * (-3 + 4 * pow(cos_2_sig, 2))))
    
    s_A_B = b * A * (sig - del_sig)   
    
    licz_ab = num.cos(U_b) * num.sin(L2)
    
    mian_ab = num.cos(U_a) * num.sin(U_b) - num.sin(U_a) * num.cos(U_b) * num.cos(L2)
    
    A_AB = num.arctan(licz_ab / mian_ab)
    
    
    if (licz_ab > 0 and mian_ab > 0):
        
        A_AB = num.rad2deg(A_AB)
        
    elif (licz_ab < 0 and mian_ab > 0):
        
        A_AB = num.rad2deg(A_AB + 2 * num.pi)
        
    elif (licz_ab > 0 and mian_ab < 0) or (licz_ab < 0 and mian_ab < 0):
        
        A_AB = num.rad2deg(A_AB + num.pi)
    
    
    licz_ba = num.cos(U_a) * num.sin(L2)
    
    mian_ba = -num.sin(U_a) * num.cos(U_b) + num.cos(U_a) * num.sin(U_b) * num.cos(L2)
    
    A_BA = num.arctan(licz_ba / mian_ba)
    
    
    if (licz_ba > 0 and mian_ba > 0):
        
        A_BA = num.rad2deg(A_BA + num.pi)
        
    elif (licz_ba < 0 and mian_ba > 0):
        
        A_BA = num.rad2deg(A_BA + 3 * num.pi)
        
    elif (licz_ba > 0 and mian_ba < 0) or (licz_ba < 0 and mian_ba < 0):
        
        A_BA = num.rad2deg(A_BA + 2 * num.pi)
    
    return (s_A_B, A_AB, A_BA)
        

def kijovi(a, e2, fi, lam, A_AB, s):
    s = s/2
    n = int(s / 1000)
    ds = 1000
    fi = num.deg2rad(fi)
    lam = num.deg2rad(lam)
    A_AB = num.deg2rad(A_AB)

    for i in range(n):
        M = (a * (1 - e2)) / pow(pow(pow((1 - e2 * num.sin(fi)), 2), 3), 0.5)
        N = a / pow((1 - (e2 * pow(num.sin(fi), 2))), 0.5)

        del_fi = ds * num.cos(A_AB) / M
        del_A_AB = num.sin(A_AB) * num.tan(fi) * ds / N

        fi_m = fi + 0.5 * del_fi
        Az_m = A_AB + 0.5 * del_A_AB

        Mi = (a * (1 - e2)) / pow(pow(pow((1 - e2 * num.sin(fi_m)), 2), 3), 0.5)
        Ni = a / pow((1 - (e2 * pow(num.sin(fi_m), 2))), 0.5)

        del_fi_m = num.cos(Az_m) * ds / Mi
        del_la_m = num.sin(Az_m) * ds / (Ni * num.cos(fi_m))
        del_Az_m = num.sin(Az_m) * num.tan(fi_m) * ds / Ni

        fi = fi + del_fi_m
        lam = lam + del_la_m
        A_AB =A_AB + del_Az_m

    ds = s % 1000
    M = (a * (1 - e2)) / pow(pow(pow(1 - e2 * num.sin(fi), 2), 3), 0.5)
    N = a / pow(1 - (e2 * pow(num.sin(fi), 2)), 0.5)

    del_fi = ds * num.cos(A_AB) / M
    del_A_AB = num.sin(A_AB) * num.tan(fi) * ds / N

    fi_m = fi + 0.5 * del_fi
    Az_m = A_AB + 0.5 * del_A_AB

    Mm = (a * (1 - e2)) / pow(pow(pow(1 - e2 * num.sin(fi_m), 2), 3), 0.5)
    Nm = a / pow(1 - (e2 * pow(num.sin(fi_m), 2)), 0.5)

    del_fi_m = num.cos(Az_m) * ds / Mm
    del_la_m = num.sin(Az_m) * ds / (Nm * num.cos(fi_m))
    del_Az_m = num.sin(Az_m) * num.tan(fi_m) * ds / Nm

    fi = fi + del_fi_m
    lam = lam + del_la_m
    A_AB = A_AB + del_Az_m

    fi = num.rad2deg(fi)
    lam = num.rad2deg(lam)
    A_AB = num.rad2deg(A_AB)
    return (fi, lam, A_AB)


def pole(fi1, lam1, fi2, lam2):
    
    e = pow(e2, 0.5)
    b = a * pow(1 - e2, 0.5)

    fi1 = num.deg2rad(fi1)
    fi2 = num.deg2rad(fi2)
    lam1 = num.deg2rad(lam1)
    lam2 = num.deg2rad(lam2)

    i1 = num.sin(fi1)/(1 - e2 * pow(num.sin(fi1), 2)) + 1/(2 * e) * num.log((1 + e * num.sin(fi1))/
                                                                              (1 - e * num.sin(fi1)))
    i2 = num.sin(fi2)/(1 - e2 * pow(num.sin(fi2), 2)) + 1/(2 * e) * num.log((1 + e * num.sin(fi2))/
                                                                              (1 - e * num.sin(fi2)))
    i = i1 - i2
    
    return (pow(b, 2) * (lam2 - lam1) * 0.5 * i)
    

x = srodek(B[0], C[0])
y = srodek(B[1], C[1])

print("Współrzędne punktu średniej szerokości:",'\n','\n', "fi:", dec_st(x),'\n', "lam:", dec_st(y))

print('\n','\n',"Pierwsze uzycie algorytmu VINCENTEGO",'\n')
s, A_ab, A_ba = vincent(a,e2,A[0], A[1], D[0], D[1])
print(" Długość linii geodezyjnej:", round(s, 5), '\n', "Azymut wprost", dec_st(A_ab), '\n', "Azymut odwrotny:", dec_st(A_ba))

print('\n','\n',"Użycie algorytmu KIJOVI")
fi, lam, Az_12 = kijovi(a, e2, A[0], A[1], A_ab, s)
print(" Współrzędne punktu środkowego:",'\n', "fi:", dec_st(fi),'\n', "lam:", dec_st(lam))

print('\n','\n', "Ponowne użycie algorytmu VINCENTEGO")
s1, A_xy, A_yx = vincent(a, e2, x, y, fi, lam)
print('\n',"Różnica odległości pomiędzy punktem środkowym a punktem średniej szerokości:", round(s1, 5),"m", '\n', "Azymut wprost:", dec_st(A_xy), '\n', "Azymut odwrotny:",
          dec_st(A_yx))

p = pole(A[0], A[1], D[0], D[1])
print('\n','\n',"Pole powierzchni:", p, "m^2")