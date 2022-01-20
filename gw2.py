import numpy as num
from tkinter import *
import matplotlib.pyplot as plt
from matplotlib import cm
from PIL import ImageTk,Image

window = Tk()
window.geometry("1080x400")
window.title("Astronomia geodezyjna ")
window.configure(background = 'black')

Warzszawa = [52.1347, 21.0042]
Nairobi = [-1.1659, 36.4900]
Sydney = [-33.5204, 151.1226]

# 8h 16min 30.924sec
rek = (( 8 + 16/60 + 30.924/3600) * 15) 
git_rek = num.deg2rad(rek)

#  9h 11min 7.98sec
dek = ( 9 + 11/60 + 7.98/3600)
git_dek = num.deg2rad(dek)

# Gwiazda z gwiazdozbioru Raka - Beta Cancri
# Wawa / Nairobi / Sydney

def greg_jul(y, m, d, h):
    if m <= 2:
        y = y - 1
        m = m + 12
        
        return (num.floor(365.25 * (y + 4716)) + num.floor(30.6001 * (m + 1)) + d + h / 24 - 1537.5)
    else:
        
        return (num.floor(365.25 * (y + 4716)) + num.floor(30.6001 * (m + 1)) + d + h / 24 - 1537.5)


def GMST(y, m, d, h):
    T = (( greg_jul(y, m, d, h) - 2451545) / 36525)
    g = ( 280.46061837 + 360.98564736629 * (greg_jul(y, m, d, h) - 2451545.0) + 0.000387933 * pow(T, 2) - pow(T, 3 / 38710000))
    
    return g%369
    


def kat_h(y, m, d, h, lam, alf):
    g = ( GMST(y, m, d, 0))
    UT1 = (h * 1.002737909350795)
    S = ( UT1 * 15 + lam + g)
    
    return ( S - alf * 15)


def zenit(fi, git_dek, t):
    fi = num.deg2rad(fi)
    t = num.deg2rad(t)
    z = num.arccos(num.sin(fi) * num.sin(git_dek) + num.cos(fi) * num.cos(git_dek) * num.cos(t))
    
    return (num.rad2deg(z))


def azymut(fi, git_dek, t):
    fi = num.deg2rad(fi)
    t = num.deg2rad(t)
    licz = -num.cos(git_dek) * num.sin(t)
    mian = num.cos(fi) * num.sin(git_dek) - num.sin(fi) * num.cos(git_dek) * num.cos(t)
    tag_A = num.rad2deg(num.arctan(licz / mian))
    
    if mian > 0 and licz < 0:
        tag_A += 360 
    elif mian < 0:
        tag_A += 180
        
    return tag_A

def wysokosc(fi, git_dek, t):
    h = 90 - zenit(fi, git_dek, t)
    
    return h

def transf_wspol(fi, git_dek, t):
    x = num.sin(num.deg2rad(zenit(fi, git_dek, t))) * num.cos(num.deg2rad(azymut(fi, git_dek, t)))
    y = num.sin(num.deg2rad(zenit(fi, git_dek, t))) * num.sin(num.deg2rad(azymut(fi, git_dek, t)))
    z = num.cos(num.deg2rad(zenit(fi, git_dek, t)))
    
    return x, y, z


godz = []

for i in range(0,24):
    godz.append(i)


godz = num.asarray(godz)


t_w = kat_h(2021, 11, 30, godz, Warzszawa[1], rek)
t_n = kat_h(2021, 11, 30, godz, Nairobi[1], rek)
t_s = kat_h(2021, 11, 30, godz, Sydney[1], rek)

# wysokosc gwiazdy nad horyzontem od czasu
h_w = wysokosc(Warzszawa[0], dek, t_w)
h_n = wysokosc(Nairobi[0], dek, t_n)
h_s = wysokosc(Sydney[0], dek, t_s)


z_w = zenit(Warzszawa[0], dek, t_w)
z_s = zenit(Sydney[0], dek, t_s)
z_n = zenit(Nairobi[0], dek, t_n)

    # odległośc zenitalna gwiazd

def w12d():
    plt.plot(godz, z_w)
    plt.xlabel('Czas (h)')
    plt.ylabel('Odległość zenitalna gwiazdy')
    plt.title('Odległość zenitalna(t) Warszawa')
    plt.show()

def w22d():
    plt.plot(godz, z_s)
    plt.xlabel('Czas (h)')
    plt.ylabel('Odległość zenitalna gwiazdy')
    plt.title('Odległość zenitalna(t) Sydney')
    plt.show()

def w32d():
    plt.plot(godz, z_n)
    plt.xlabel('Czas (h)')
    plt.ylabel('Odległość zenitalna gwiazdy')
    plt.title('Odległość zenitalna(t) Nairobi')
    plt.show()

    #wykresy wysokości gwiazdy nad horyzontem
def w42d():
    plt.plot(godz, h_w)
    plt.xlabel('Czas (h)')
    plt.ylabel('Wysokośc gwiazdy nad horyzontem')
    plt.title('H(t) Warszawa')
    plt.show()

def w52d():
    plt.plot(godz, h_s)
    plt.xlabel('Czas (h)')
    plt.ylabel('Wysokośc gwiazdy nad horyzontem')
    plt.title('H(t) Sydney')
    plt.show()

def w62d():
    plt.plot(godz, h_n)
    plt.xlabel('Czas (h)')
    plt.ylabel('Wysokośc gwiazdy nad horyzontem')
    plt.title('H(t) Nairobi')
    plt.show()


b = num.vectorize(azymut)
c = num.vectorize(wysokosc)
d = num.vectorize(transf_wspol)
x_wa, y_wa, z_wa = d(Warzszawa[0], git_dek, t_w)
x_pm, y_pm, z_pm = d(Sydney[0], git_dek, t_s)
x_n, y_n, z_n = d(Nairobi[0], git_dek, t_n)

u, v = num.mgrid[0:2 * num.pi:30j, 0:num.pi:20j]
x = num.cos(u) * num.sin(v)
y = num.sin(u) * num.sin(v)
z = num.cos(v)

def w13d():
    # fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_title("Ruch gwiazdy po niebie (Warszawa)")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.plot_surface(x, y, z, alpha=0.05, facecolors=cm.jet(z/num.amax(z)))
    ax.plot3D(x_wa, y_wa, z_wa, 'blue')
    ax.scatter3D(x_wa, y_wa, z_wa, c='none', cmap='cividis')
    plt.show()

def w23d():
    ax = plt.axes(projection='3d')
    ax.set_title("Ruch gwiazdy po niebie (Nairobi)")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.plot_surface(x, y, z, alpha=0.05, facecolors=cm.jet(z / num.amax(z)))
    ax.plot3D(x_n, y_n, z_n, 'blue')
    ax.scatter3D(x_n, y_n, z_n, c='none', cmap='cividis')
    plt.show()

def w33d():
    ax = plt.axes(projection='3d')
    ax.set_title("Ruch gwiazdy po niebie (Sydney)")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.plot_surface(x, y, z, alpha=0.05, facecolors=cm.jet(z / num.amax(z)))
    ax.plot3D(x_pm, y_pm, z_pm, 'blue')
    ax.scatter3D(x_pm, y_pm, z_pm, c='none', cmap='cividis')
    plt.show()

obl = Button(window,
            text="Tor ruchu gwiazdy nad Sydney",
            command=w33d,
            font=("Arial",10,'bold'),
            activeforeground="black",
            activebackground="#F0FFC0")
obl.place(x=820, y=250)

obl = Button(window,
            text="Tor ruchu gwiazdy nad Nairobi",
            command=w23d,
            font=("Arial",10,'bold'),
            activeforeground="black",
            activebackground="#F0FFC0")
obl.place(x=450, y=250)

obl = Button(window,
            text="Tor ruchu gwiazdy nad Warszawą",
            command=w13d,
            font=("Arial",10,'bold'),
            activeforeground="black",
            activebackground="#F0FFC0")
obl.place(x=100, y=250)

obl = Button(window,
            text="Odległość zenitalna gwiazdy od Warszawy",
            command=w12d,
            font=("Arial",10,'bold'),
            activeforeground="black",
            activebackground="#F0FFC0")
obl.place(x=70, y=80)

obl = Button(window,
            text="Odległość zenitalna gwiazdy od Warszawy Sydney",
            command=w22d,
            font=("Arial",10,'bold'),
            activeforeground="black",
            activebackground="#F0FFC0")
obl.place(x=390, y=80)

obl = Button(window,
            text="Odległość zenitalna gwiazdy od Nairobi",
            command=w32d,
            font=("Arial",10,'bold'),
            activeforeground="black",
            activebackground="#F0FFC0")
obl.place(x=790, y=80)

obl = Button(window,
            text="Wysokośc gwiazdy nad horyzontem (Warszawa)",
            command=w42d,
            font=("Arial",10,'bold'),
            activeforeground="black",
            activebackground="#F0FFC0")
obl.place(x=40, y=170)

obl = Button(window,
            text="Wysokośc gwiazdy nad horyzontem (Sydney)",
            command=w52d,
            font=("Arial",10,'bold'),
            activeforeground="black",
            activebackground="#F0FFC0")
obl.place(x=430, y=170)

obl = Button(window,
            text="Wysokośc gwiazdy nad horyzontem (Nairobi)",
            command=w62d,
            font=("Arial",10,'bold'),
            activeforeground="black",
            activebackground="#F0FFC0")
obl.place(x=780, y=170)


window.mainloop()