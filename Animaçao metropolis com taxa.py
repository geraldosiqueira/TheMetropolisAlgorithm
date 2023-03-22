# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 09:21:02 2022

@author: Geraldo Siqueira
"""

import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.animation import FuncAnimation

k = 1
Kb = 1

def U(x,y):
    return (1/2)*k*(x**2 + y**2)

def Histograma(x, t, Nhist, a, b):
    x_aux = []
    for k in range( np.size(x[:,t])):
        if x[k,t] >= a and x[k,t] <= b:
            x_aux.append(x[k,t])
        
    dx = (b-a)/Nhist
    Xh = np.arange(a, b, dx)
    Hist = np.zeros(Nhist)
    N = np.size(x_aux)
    for i in range(N):
        j = int((x_aux[i] - a)/dx)
        if j<Nhist and j>0:
            Hist[j] += 1
    Hist = Hist/(N*dx)
    return Xh, Hist

T = 0.1
Nt = 1000
N = 1000
D = 0.7
k = 1

L = 2
A = np.sqrt(k/(2*np.pi*Kb*T))
X = np.arange(-L, L, 0.01)
P = A*np.exp(-k*(X**2)/(2*Kb*T))

x = np.zeros( [N, Nt] )
y = np.zeros( [N, Nt] )
taxa = np.zeros(Nt)

#Fazendo as tentativas de mudanças nos valores de x e y:
for n in range(0,Nt-1):
    vp_x = np.zeros(N)
    vp_y = np.zeros(N)
    aux = 0
    for i in range (0,N-1):
        Part = random.randint(0,N-1)
        varx = random.uniform(-D,D)
        vary = random.uniform(-D,D)
        x[:,n+1] = x[:,n]
        y[:,n+1] = y[:,n]
        x[Part, n+1] += varx
        y[Part, n+1] += vary
        
        DeltaE = U(x[Part, n+1] + vp_x[Part], y[Part, n+1] + vp_y[Part]) - U(x[Part, n] + vp_x[Part], y[Part, n] + vp_y[Part])
        x[Part,n+1] = x[Part,n]
        y[Part,n+1] = y[Part,n]
        #print("ΔE=", DeltaE)
        if DeltaE > 0:
            u = random.uniform(0,1)
            p = np.exp((-DeltaE)/(Kb*T))
            if u > p:
                varx, vary = 0, 0
                aux += 1
        vp_x[Part] += varx
        vp_y[Part] += vary
        
    x[:,n+1] += vp_x[:]
    y[:,n+1] += vp_y[:]
    taxa[n] = (N - aux)/N
    
dt = 0.01

#Configurando a figura onde será feita a animação
L = 1.0
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(autoscale_on=False, xlim=(-L, L), ylim=(-L, L))
ax.set_aspect('equal')
ax.set_title("Metropole com Aceitação (Δ = 0.7)")

# Objetos que receberão os dados (configurações) e texto a serem mostrados em cada quadro
confs, = ax.plot([], [], 'ko', color = 'black', ms=5, alpha=0.2)
texto = ax.text(0.05, 0.9, '', transform=ax.transAxes)

# Converte os dados recebidos em cada iteração em objetos a serem mostrados na figura
def animate(n):
    confs.set_data(x[:,n], y[:,n])
    texto.set_text('Taxa de Aceitação = %.2f' % taxa[n])
    return confs, texto

# Constrói a animação
ani = FuncAnimation(fig, animate, frames=Nt,interval=50)


plt.show()
ani.save("metropole_com_taxa_0_7.gif")