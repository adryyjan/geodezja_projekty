from errno import EMSGSIZE
import numpy as np


global kappa
global a 
global e2
global x0 
global y0 
global z0 
global ex 
global ey 
global ez
global h  
 

a = 6378137 
e2 = 0.0066943800290

a1 = 6378245
e21 = 0.0066934215520

x0 = -33.4297
y0 = 146.5746
z0 = 76.2865

ex = -0.35867/3600
ey = -0.05283/3600
ez = 0.84354/3600

kappa = pow(0.8407728 * 10, -6)

A = [50.25, 20.75]
B = [50, 20.75]
C = [50.25, 21.25]
D = [50, 21.25]

h = 1

def dec_st(x):
    
    h = int(x)
    m = int((x - h) * 60)
    s = round((x - h - m/60) * 3600, 5)

    return (str(str(h) + "°"+str(m) + "'"+str(s) + "''"))

def zmiana( s, m = 0.0 ,sec = 0.0):
    
    r = s + m/60 + sec/3600
    
    return r


def do_xyz(fi, lam, h):
    
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = a / np.sqrt(1 - e2 * np.sin(fi) ** 2)
    x = (N + h) * np.cos(fi) * np.cos(lam)
    y = (N + h) * np.cos(fi) * np.sin(lam)
    z = (N * (1 - e2) + h) * np.sin(fi)
    return x, y, z

def hirv(x, y, z, a, e2):
    
    r = (x ** 2 + y ** 2) ** 0.5
    fi = np.arctan((z/r) * pow((1-e2), -1))

    N = a/np.sqrt(1-e2 * pow(np.sin(fi), 2))
    h = r/np.cos(fi) - N
    fi2 = np.arctan((z / r) * pow((1 - e2 * (N / (N + h))), -1))
    sek = np.deg2rad(0.00005/3600)
    
    while abs(fi2-fi) >= sek:
        
        fi = fi2
        N = a / np.sqrt(1 - e2 * pow(np.sin(fi), 2))
        h = r / np.cos(fi) - N
        fi2 = np.arctan((z / r) * pow((1 - e2 * (N / (N + h))), -1))
    N = a / np.sqrt(1 - e2 * pow(np.sin(fi2), 2))
    h = r / np.cos(fi2) - N
    lam = np.arctan(y/x)

    fi2 = np.rad2deg(fi2)
    lam = np.rad2deg(lam)

    return fi2, lam, h

def transform(x, y, z):
    
    al = np.deg2rad(ex)
    bet = np.deg2rad(ey)
    gm = np.deg2rad(ez)
    
    M = np.array([
        [kappa, gm, -bet],
        [-gm, kappa, al],
        [bet, -al, kappa]])
    
    M_p = np.array([x,y,z])
    
    M_0 = np.array([x0,y0,z0])
    
    return (M_p + M @ M_p + M_0)


srodkowy = []

srodkowy.append(zmiana(49,59,33.94))
srodkowy.append(zmiana(21,0,6.62))


sr_szer = [] 

sr_szer.append(zmiana(50,7,30))
sr_szer.append(zmiana(21))

xa, ya, za = do_xyz(A[0], A[1], h)
xb, yb, zb = do_xyz(B[0], B[1], h)
xc, yc, zc = do_xyz(C[0], C[1], h)
xd, yd, zd = do_xyz(D[0], D[1], h)
x_srodkowy, y_srodkowy, z_srodkowy = do_xyz(srodkowy[0], srodkowy[1], h)
xsr_szer, ysr_szer, zsr_szer = do_xyz(sr_szer[0], sr_szer[1], h)

fi_A, lam_A, h_A = transform(xa, ya, za)
fi_B, lam_B, h_B = transform(xb, yb, zb)
fi_C, lam_C, h_C = transform(xc, yc, zc)
fi_D, lam_D, h_D = transform(xd, yd, zd)
fi_srodkowy, lam_srodkowy, h_srodkowy = transform(x_srodkowy, y_srodkowy, z_srodkowy)
fi_sr_szero, lam_sr_szero, h_sr_szero = transform(xsr_szer, ysr_szer, zsr_szer)

fi_a, lam_a, h_a = hirv(fi_A, lam_A, h_A, a1, e21)
fi_b, lam_b, h_b = hirv(fi_B, lam_B, h_B, a1, e21)
fi_c, lam_c, h_c = hirv(fi_C, lam_C, h_C, a1, e21)
fi_d, lam_d, h_d = hirv(fi_D, lam_D, h_D, a1, e21)
fi_srod, lam_srod, h_srod = hirv(fi_srodkowy, lam_srodkowy, h_srodkowy, a1, e21)
fi_sr_szer, lam_sr_szer, h_sr_szer = hirv(fi_sr_szero, lam_sr_szero, h_sr_szero, a1, e21)

    #zamiana dzies na stopnie
fi_a = dec_st(float(fi_a))
fi_b = dec_st(float(fi_b))
fi_c = dec_st(float(fi_c))
fi_d = dec_st(float(fi_d))
lam_a = dec_st(float(lam_a))
lam_b = dec_st(float(lam_b))
lam_c = dec_st(float(lam_c))
lam_d = dec_st(float(lam_d))

fi_srodkowa = dec_st(float(fi_srod))
fi_sr_szer = dec_st(float(fi_sr_szer))
lam_sr_szer = dec_st(float(lam_sr_szer))
lam_srodkowa = dec_st(float(lam_srod))

print(
    "A:", transform(xa, ya, za), '\n',
    "B:", transform(xb, yb, zb), '\n',
    "C:", transform(xc, yc, zc), '\n',
    "D:", transform(xd, yd, zd), '\n',
    "Środkowy:", transform(fi_srodkowy, lam_srodkowy, h_srodkowy), '\n',
    "Średniej szerokości:", transform(fi_sr_szero, lam_sr_szero, h_sr_szero))

print('\n', "_____________________________________________", '\n')
print(
    " Punkt A:", fi_a, lam_a, round(h_a, 6), '\n',
    "Punkt B:", fi_b, lam_b, round(h_b, 6), '\n',
    "Punkt C:", fi_c, lam_c, round(h_c, 6), '\n',
    "Punkt D:", fi_d, lam_d, round(h_d, 6), '\n',
    "Środkowy:", fi_sr_szer, lam_srodkowy,round(h_srodkowy,6), '\n',
    "Średniej szerokości:", fi_sr_szer, lam_sr_szer, round(h_sr_szer, 6), '\n')


    