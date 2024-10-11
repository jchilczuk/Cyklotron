import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle


def F_net(r, v, t, B, q, m, V):
    d = 90e-6  # przerwa między duantaki
    E_0 = V / (d)  # pole elektryczne między duantami

    w = q * np.linalg.norm(B) / m  # czestotliwość cyklotronu

    F = np.zeros(3)
    if np.absolute(r[0]) < d / 2:
        F[0] = q * E_0 * np.cos(w * t)
    else:
        F = q * np.cross(v, B)
    return F


def cyclotron(q, m, B, dt, T, V):
    N = int(T / dt)
    r_cyclotron = .05  # promień duantów

    r = np.zeros((N, 3))
    v = np.zeros((N, 3))
    t = np.zeros(N)

    F = np.zeros((N, 3))

    i = 0
    while (np.linalg.norm(r[i]) < r_cyclotron):
        F[i] = F_net(r[i], v[i], t[i], B, q, m, V)
        v[i + 1] = v[i] + F[i] / m * dt
        r[i + 1] = r[i] + v[i + 1] * dt
        t[i + 1] = t[i] + dt
        if abs(np.linalg.norm(v[i+1])) > 3e8:
            print("Cząstka osiągnęła prędkość większą od prędkości światła.")
            sys.exit(0)
        else:
            i += 1

    Ek = 0.5*m*(np.linalg.norm(v[i]))**2
    print(f"Prędkość cząstki po opuszczeniu cyklotronu: {np.linalg.norm(v[i]):.3e} m/s ={np.linalg.norm(v[i])/3e8: .3f} prędkości światła")
    print(f"Energia kinetyczna cząstki po opuszczeniu cyklotronu: {Ek/q:.3e} eV")
    print(f"Całkowity czas cząstki w cyklotronie: {t[i]:.3e} s" )

    return r[:i], v[:i], t[:i]


print("\nWitaj, przedstawiamy Ci symulację przyspieszania cząstek w cyklotronie.\n")
print("Wybierz rodzaj cząstki: 1 - proton, 2 - cząstka alfa, 3 - deuteron")
x = input("Numer wybranej cząstki: ")

if x == '1':
    q = 1.6e-19
    m = 1.67e-27
elif x == '2':
    q = 2*1.6e-19
    m = 4*1.67e-27
elif x == '3':
    q = 1.6e-19
    m = 2*1.67e-27
else:
    print("Nie ma takiej cząstki")
    sys.exit(0)

b = float(input("Podaj wartość indukcji magnetycznej B [T]: "))
V = float(input("Podaj wartość napięcia między duantami [V]: "))
print("\nZaczynamy!\n")
B = np.array([0.0, 0.0, b])

dt = 5e-12
T = 1e-5
d = 90e-6
r, v, t = cyclotron(q, m, B, dt, T, V)

fig, ax = plt.subplots()
ax.set_xlim(-0.05, 0.05)
ax.set_ylim(-0.05, 0.05)
ax.set_aspect('equal')
circle = Circle((0, 0), 0.05, edgecolor='c', facecolor='c')
ax.add_patch(circle)
ax.axvline(x=0, color='w')
fig.set_size_inches(8,7)
fig.canvas.manager.window.move(0,50)
line, = ax.plot([], [], 'b', linestyle='--' )

def init():
    line.set_data([], [])
    return line,


def update(frame):
    step = len(r) // 100
    line.set_data(r[:frame*step, 0], r[:frame*step, 1])
    return line,

ani = FuncAnimation(fig, update, frames=100, init_func=init, blit=True)

plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('Ruch cząstki w cyklotronie')

#plt.get_current_fig_manager().window.geometry("+0+50")
#fig.set_size_inches(8,6)
# plt.show()

fig2, ax2 = plt.subplots()
ax2.set_xlim(0,t[-1])
ax2.set_ylim(0,np.linalg.norm(v, axis=1)[-1])
fig2.set_size_inches(9,7)
fig2.canvas.manager.window.move(800,50)
line2, = ax2.plot([], [], 'b-')

def init():
    line2.set_data([], [])
    return line2,

v2 = np.linalg.norm(v, axis=1)

def update(frame):
    step = len(v2)// 100
    line2.set_data(t[:frame*step], v2[:frame*step])
    return line2,

ani2 = FuncAnimation(fig2, update, frames=100, init_func=init, blit=True)

plt.xlabel('t [s]')
plt.ylabel('v [m/s]')
plt.title('Prędkość cząstki w cyklotronie')
#fig2.set_size_inches(8, 6)
#plt.get_current_fig_manager().window.geometry("+800+50")

plt.show()



