import numpy as num
from shapely.geometry import Polygon


global a 
global e2
global m_20 
global m_92 

a = 6378137 
e2 = 0.0066943800290 

m_92 = 0.9993
m_20 = 0.999923

A = [50.25, 20.75]
B = [50, 20.75]
C = [50.25, 21.25]
D = [50, 21.25]



def zmiana( s, m = 0.0 ,sec = 0.0):
    
    r = s + m/60 + sec/3600
    
    return r


def do_GK(fi, lam, lam0):
    fi = num.deg2rad(fi)
    lam = num.deg2rad(lam)

    b2 = a ** 2 * (1 - e2)
    e2_prim = (pow(a, 2) - b2) / b2
    lam0 = num.deg2rad(lam0)
    t = num.tan(fi)
    N = a / pow(1 - e2 * pow(num.sin(fi), 2),0.5)
    eta2 = e2_prim * pow(num.cos(fi), 2)
    delta_lam = lam - lam0

    A0 = 1 - e2 / 4 - 3 * pow(e2, (2 / 64)) - 5 * pow(e2, (3 / 256))
    A2 = 3 / 8 * (e2 + pow(e2, (2 / 4)) + 15 * pow(e2, (3 / 128)))
    A4 = 15 / 256 * (pow(e2, 2) + 3 * pow(e2, (3 / 4)))
    A6 = 35 * pow(e2, (3 / 3072))

    sig = a * (A0 * fi - A2 * num.sin(2 * fi) + A4 * num.sin(4 * fi) - A6 * num.sin(6 * fi))

    x = sig + pow(delta_lam, (2 / 2)) * N * num.sin(fi) * num.cos(fi) * (1 + pow(delta_lam, (2 / 12)) * pow(num.cos(fi), 2) * (5 -
        pow(t, 2) + 9 * eta2 + 4 * pow(eta2, 2)) + pow(delta_lam, (4 / 360)) * pow(num.cos(fi), 4) * (61 - 58 * pow(t, 2) + pow(t, 4) +
        270 * eta2 - 330 * eta2 * pow(t, 2)))

    y = delta_lam * N * num.cos(fi) * (
                1 + pow(delta_lam, (2 / 6)) * pow(num.cos(fi), 2) * (1 - pow(t, 2) + eta2) + pow(delta_lam, (4 / 120))
                * pow(num.cos(fi), 4) * (5 - 18 * pow(t, 2) + pow(t, 4) + 14 * eta2 - 58 * eta2 * pow(t, 2)))
    return x, y


def z_GK(x, y, lam0):
    lam0 = num.deg2rad(lam0)
    a2 = a ** 2
    b2 = a2 * (1 - e2)
    e2_prim = (a2 - b2) / b2

    A0 = 1 - e2 / 4 - 3 *pow(e2 , (2 / 64)) - 5 * pow(e2, (3 / 256))
    A2 = 3 / 8 * (e2 + pow(e2, (2 / 4)) + 15 * e2 ** 3 / 128) #dlacego pow tutaj wywala prgoram????
    A4 = 15 / 256 * (pow(e2, 2) + 3 * pow(e2, (3 / 4)))
    A6 = 35 * e2 ** 3 / 3072

    fi = x / (a * A0)
    sig = a * (A0 * fi - A2 * num.sin(2 * fi) + A4 * num.sin(4 * fi) - A6 * num.sin(6 * fi))

    while True:
        fi2 = fi + (x - sig) / (a * A0)

        sig = a * (A0 * fi2 - A2 * num.sin(2 * fi2) + A4 * num.sin(4 * fi2) - A6 * num.sin(6 * fi2))
        N = a / pow((1 - e2 * pow(num.sin(fi2), 2)), 0.5)
        M = (a * (1 - e2)) / pow((pow(1 - e2 * pow(num.sin(fi2), 2), 3)), 0.5)
        t = num.tan(fi2)
        eta2 = e2_prim * pow(num.cos(fi2), 2)

        if abs(fi2 - fi) < num.deg2rad(0.000001 / 3600):
            break

        fi = fi2

    fi = fi2 - pow(y, 2) * t / (2 * M * N) * (1 - y ** 2 / (12 * pow(N, 2)) * (5 + 3 * pow(t, 2) + eta2 - 9 * eta2 *
        pow(t, 2) - 4 * pow(eta2, 2)) + y ** 4 / (360 * pow(N, 4)) * (61 + 90 * pow(t, 2) + 45 * pow(t, 4))) #tuataj tez pow( nie działa na długim wierszu)

    lam = lam0 + y / (N * num.cos(fi2) * (1 - y ** 2 / (6 * pow(N, 2)) * (1 + 2 * pow(t, 2) + eta2) + y ** 4 / (120 * pow(N, 4))
        * (5 + 28 * pow(t, 2) + 24 * pow(t, 4) + 6 * eta2 + 8 * eta2 * pow(t, 2)))) #tutaj również :(((
            
    fi = num.rad2deg(fi)
    lam = num.rad2deg(lam)
    return fi, lam


def do_1992(x, y):
    
    x = m_92 * x - 5300000
    y = m_92 * y + 500000
    return x, y

def z_1992(x, y):
    
    x = (x + 5300000)/m_92
    y = (y - 500000)/m_92
    return x, y


def do_2000(x, y, lam):
    lam = num.deg2rad(lam)
    
    
    if num.deg2rad(22.5) <= lam <= num.deg2rad(25.5) :
        nr = 8
    elif num.deg2rad(19.5) <= lam <= num.deg2rad(22.5) :
        nr = 7
    elif num.deg2rad(16.5) <= lam <= num.deg2rad(19.5):
        nr = 6
    elif num.deg2rad(13.5) <= lam <= num.deg2rad(16.5):
        nr = 5

    x = m_20 * x
    y = m_20 * y + nr * 1000000 + 500000

    return x, y

def z_2000(x, y, nr):
    
    x = x/m_20
    y = (y - 500000 - nr * 1000000)/m_20
    return x, y

def skala_1992(x_gk, y_gk):
    x, y = z_1992(x_gk, y_gk)
    fi= z_GK(x, y, 19)
    M = a * (1 - e2)/(1 - e2 * pow(pow(num.sin(fi), 2), (3/2)))
    N = a / pow((1 - e2 * pow(num.sin(fi), 2)),  0.5)

    Q = num.sqrt(M * N)

    m_gk = 1 + y ** 2/(2 * pow(Q, 2)) + y ** 2/(24 * pow(Q, 4)) # też...
    
    m_92 = 0.9993 * m_gk
    kappa = (1 - m_92) * 1000

    return m_92, kappa

def skala_2000(x_gk, y_gk):
    
    x, y = z_2000(x_gk, y_gk, 7)
    fi= z_GK(x, y, 21)
    M = a * (1 - e2) / pow((1 - e2 * pow(num.sin(fi), 2)), (3 / 2))
    N = a / pow((1 - e2 * pow(num.sin(fi), 2)), 0.5)
    Q = num.sqrt(M * N)

    m_gk = 1 + y ** 2 / (2 * pow(Q, 2)) + y ** 2 / (24 * pow(Q, 4))
    m2000 = m_20 * m_gk
    kappa = (1 - m2000) * 1000

    return m2000, kappa

def skala_GK(x_gk, y_gk):
    fi = z_GK(x_gk, y_gk, 19)
    M = a * (1 - e2) / pow((1 - e2 * pow(num.sin(fi), 2)), (3 / 2))
    N = a / pow((1 - e2 * pow(num.sin(fi), 2)), 0.5)
    Q = num.sqrt(M * N)
    m_gk = 1 + y_gk ** 2 / (2 * pow(Q, 2)) + y_gk ** 2 / (24 * pow(Q, 4))
    kappa = (1 - m_gk) * 1000
    return m_gk, kappa

def skala_pola(m, kappa):
    m = pow(m, 2)
    kappa = (1 - m) * 10000
    return m, kappa


srodkowy = []

srodkowy.append(zmiana(49,59,33.94))
srodkowy.append(zmiana(21,0,6.62))


sr_szer = [] 

sr_szer.append(zmiana(50,7,30))
sr_szer.append(zmiana(21))

x_gk = []
y_gk = []

x_1992 = []
y_1992 = []

x_2000 = []
y_2000 = []

tab_do_pol = num.array([])
tab_do_pol_1 = [A[0],B[0],C[0],D[0],A[0],sr_szer[0],srodkowy[0]]

tab_do_pol = num.array([])
tab_do_pol_2 = [A[1],B[1],C[1],D[1],A[1],sr_szer[1],srodkowy[1]]


for i in range(0,7):
    x1_kg, y1_kg = do_GK(tab_do_pol_1[i], tab_do_pol_2[i], 19)
    x_gk.append(x1_kg)
    y_gk.append(y1_kg)
    x1_1992, y1_1992 = do_1992(x1_kg, y1_kg)
    x_1992.append(x1_1992)
    y_1992.append(y1_1992)
    x1_2000, y1_2000 = do_2000(x1_kg, y1_kg, D[1])
    x_2000.append(x1_2000)
    y_2000.append(y1_2000)


px_gk = num.array([])
px_gk = [x_gk[0],x_gk[1],x_gk[3],x_gk[2],x_gk[1]]
py_gk = num.array([])
py_gk = [y_gk[0],y_gk[1],y_gk[3],y_gk[2],y_gk[1]]



wsp_2000 = ((5569691.697548566, 7624802.61719599), (5541888.267106318, 7625454.540104091), (5542847.059764723, 7661295.862670981), (5570648.987065094, 7660457.55048384), (5569691.697548566, 7624802.61719599))
p_20 = Polygon(wsp_2000)
pole_2000 = p_20.area/1000000

wsp_1992 = ((266221.5124167381, 624724.8591781091), (238435.4048455162, 625376.375906963), (239393.60013009794, 661195.367610418), (267178.20549597126, 660357.5777319864), (266221.5124167381, 624724.8591781091))
p_92 = Polygon(wsp_1992)
pole_1992 = p_92.area/1000000

wsp_gk = ((5570120.596834523, 124812.22773752533), (5542315.025363271, 125464.20084755628), (5543273.891854396, 161308.28340880413),(5571077.960068019, 160469.9066666531), (5570120.596834523, 124812.22773752533))
p_gk = Polygon(wsp_gk)
pole_gk = p_gk.area/1000000


mGK_A, kappaGK_A = skala_GK(x_gk[0], y_gk[0])
mGK_B, kappaGK_B = skala_GK(x_gk[1],y_gk[1])
mGK_C, kappaGK_C = skala_GK(x_gk[2], y_gk[2])
mGK_D, kappaGK_D = skala_GK(x_gk[3], y_gk[3])
mGK_srszer, kappaGK_srszer = skala_GK(x_gk[4], y_gk[4])
mGK_srodkowy, kappaGK_srodkowy = skala_GK(x_gk[5],y_gk[5])


m1992_A, kappa1992_A = skala_1992(x_1992[0], y_1992[0])
m1992_B, kappa1992_B = skala_1992(x_1992[1], y_1992[1])
m1992_C, kappa1992_C = skala_1992(x_1992[2], y_1992[2])
m1992_D, kappa1992_D = skala_1992(x_1992[3], y_1992[3])
m1992_srodkowy, kappa1992_srodkowy = skala_1992(x_1992[4], y_1992[4])
m1992_srszer, kappa1992_srszer = skala_1992(x_1992[4], y_1992[4])

m2000_A, kappa2000_A = skala_2000(x_2000[0], y_2000[0])
m2000_B, kappa2000_B = skala_2000(x_2000[1], y_2000[1])
m2000_C, kappa2000_C = skala_2000(x_2000[2], y_2000[1])
m2000_D, kappa2000_D = skala_2000(x_2000[3], y_2000[3])
m2000_srszer, kappa2000_srszer = skala_2000(x_2000[4],y_2000[4])
m2000_srodkowy, kappa2000_srodkowy = skala_2000(x_2000[4],y_2000[4])


GK_m_A, GK_kappa_A = skala_pola(mGK_A, kappaGK_A)
GK_m_B, GK_kappa_B= skala_pola(mGK_B, kappaGK_B)
GK_m_C, GK_kappa_C = skala_pola(mGK_C, kappaGK_C)
GK_m_D, GK_kappa_D= skala_pola(mGK_D, kappaGK_D)
GK_m_srszer, GK_kappa_srszer = skala_pola(mGK_srszer, kappaGK_srszer)
GK_m_srodkowy, GK_kappa_srodkowy = skala_pola(mGK_srodkowy, kappaGK_srodkowy)


m_A_1992, kappa_A_1992 = skala_pola(m1992_A, kappa1992_A)
m_B_1992, kappa_B_1992 = skala_pola(m1992_B, kappa1992_B)
m_C_1992, kappa_C_1992 = skala_pola(m1992_C, kappa1992_C)
m_D_1992, kappa_D_1992 = skala_pola(m1992_D, kappa1992_D)
m_srszer_1992, kappa_srszer_1992 = skala_pola(m1992_srszer, kappa1992_srszer)
m_srodkowy_1992, kappa_srodkowy_1992 = skala_pola(m1992_srodkowy, kappa1992_srodkowy)

m_A_2000, kappa_A_2000 = skala_pola(m2000_A, kappa2000_A)
m_B_2000, kappa_B_2000 = skala_pola(m2000_B, kappa2000_B)
m_C_2000, kappa_C_2000 = skala_pola(m2000_C, kappa2000_C)
m_D_2000, kappa_D_2000 = skala_pola(m2000_D, kappa2000_D)
m_srszer_2000, kappa_srszer_2000 = skala_pola(m2000_srszer, kappa2000_srszer)
m_srodkowy_2000, kappa_srodkowy_2000 = skala_pola(m2000_srodkowy, kappa2000_srodkowy)


print("---------------------tabela-1--------------------")
print("-----------------------GK------------------------")
print("A: " + str(x_gk[0]) + " " + str(y_gk[0]))
print("B: " + str(x_gk[1]) + " " + str(y_gk[1]))
print("C: " + str(x_gk[2]) + " " + str(y_gk[2]))
print("D: " + str(x_gk[3]) + " " + str(y_gk[3]))
print("S: " + str(x_gk[4]) + " " + str(y_gk[4]))
print("M: " + str(x_gk[5]) + " " + str(y_gk[5]))
print("-----------------------92------------------------")
print("A: " + str(x_1992[0]) + " " + str(y_1992[0]))
print("B: " + str(x_1992[1]) + " " + str(y_1992[1]))
print("C: " + str(x_1992[2]) + " " + str(y_1992[2]))
print("D: " + str(x_1992[3]) + " " + str(y_1992[3]))
print("S: " + str(x_1992[4]) + " " + str(y_1992[4]))
print("M: " + str(x_1992[5]) + " " + str(y_1992[5]))
print("----------------------2000-----------------------")
print("A: " + str(x_2000[0]) + " " + str(y_2000[0]))
print("B: " + str(x_2000[1]) + " " + str(y_2000[1]))
print("C: " + str(x_2000[2]) + " " + str(y_2000[2]))
print("D: " + str(x_2000[3]) + " " + str(y_2000[3]))
print("S: " + str(x_2000[4]) + " " + str(y_2000[4]))
print("M: " + str(x_2000[5]) + " " + str(y_2000[5]))
print("---------------------tabela-2--------------------")
print("Pole Elipsoidy: " + "994.265196074311" + " km^2.")
print("Pole GK:        " + str(pole_gk ) + " km^2.")
print("Pole 1992:      " + str(pole_1992 ) + " km^2.")
print("Pole 2000:      " + str(pole_2000 ) + " km^2.")
print("--------------------tabela-3-i-4-----------------")
print("-----------------------GK------------------------")
print("A: " + str(GK_m_A))
print("B: " + str(GK_m_B))
print("C: " + str(GK_m_C))
print("D: " + str(GK_m_D))
print("S: " + str(GK_m_srszer))
print("M: " + str(GK_m_srodkowy))
print("-----------------------92------------------------")
print("A: " + str(m_A_1992))
print("B: " + str(m_B_1992))
print("C: " + str(m_C_1992))
print("D: " + str(m_D_1992))
print("S: " + str(m_srszer_1992))
print("M: " + str(m_srodkowy_1992))
print("----------------------2000-----------------------")
print("A: " + str(m_A_2000))
print("B: " + str(m_B_2000))
print("C: " + str(m_C_2000))
print("D: " + str(m_D_2000))
print("S: " + str(m_srszer_1992))
print("M: " + str(m_srodkowy_1992))
print("-------------------------------------------------")



