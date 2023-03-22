# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 14:46:43 2022

@author: Geraldo Siqueira
"""

import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.animation import FuncAnimation
from numba import jit

k = 1
Kb = 1
a = 1
e = 1

@jit
def U(xi, yi, xj, yj):
    r2 = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj)
    r6 = r2*r2*r2
    r12 = r6*r6
    return 4*e*((1/r12) - (1/r6))


Nt = 600000
N = 36
D = 0.002*a
k = 1

L = 2
'''
A = np.sqrt(k/(2*np.pi*Kb*T))
X = np.arange(-L, L, 0.01)
P = A*np.exp(-k*(X**2)/(2*Kb*T))
'''

x = np.zeros( [N, Nt] )
y = np.zeros( [N, Nt] )

x[0,0] = random.uniform(-L,L)
y[0,0] = random.uniform(-L,L)
for l in range(1,N-1):
    for g in range(N):
        if l != g:
            while ((x[l,0]-x[g,0])*(x[l,0]-x[g,0]) + (y[l,0]-x[g,0])*(y[l,0]-x[g,0]))**2 < a :
                x[l,0] = random.uniform(-L,L)
                y[l,0] = random.uniform(-L,L)

#Fazendo as tentativas de mudanças nos valores de x e y:
for n in range(0,Nt-1):
    energiaF = 0
    energiaP = 0
    Part = random.randint(0,N-1)
    varx = random.uniform(-D,D)
    vary = random.uniform(-D,D) 
    x[:,n+1] = x[:,n]
    y[:,n+1] = y[:,n]
    x[Part, n+1] += varx
    y[Part, n+1] += vary
    for k in range (N):
        if Part != k:
            energiaF += U(x[Part,n+1], y[Part, n+1], x[k,n+1], y[k,n+1])
            energiaP += U(x[Part,n], y[Part, n], x[k,n], y[k,n])
    DeltaE = (1/2)*(energiaF - energiaP)
    x[Part,n+1] = x[Part,n]
    y[Part,n+1] = y[Part,n]
    #print("ΔE=", DeltaE)
    if DeltaE > 0:
        #u = random.uniform(0,1)
        #p = np.exp((-DeltaE)/(Kb*T))
        #if u > p:
        varx, vary = 0, 0
        
    x[Part,n+1] += varx
    y[Part,n+1] += vary
'''
for k in range(0,N-1):
    f = open("posicaofinalX.txt", "a")
    g = open("posicaofinalY.txt", 'a')
    f.write(str(x[k, Nt-1]) + ";")   
    g.write(str(y[k, Nt-1]) + ";")
    f.close()
    g.close()
'''

x2 = np.zeros( [N, Nt] )
y2 = np.zeros( [N, Nt] )
''' 
f = open('posicaofinalX.txt','r')
g = open('posicaofinalY.txt','r')
pfX = f.read()
pfY = g.read()
pfX_array = pfX.split(';')
pfY_array = pfY.split(';')
for m in range(0, N-1):
    pfX_array[m] = float(pfX_array[m])
    pfY_array[m] = float(pfY_array[m])
#print(pfX_array)
#print(pfY_array)

for j in range (0,N):
    x2[j, 0] = pfX_array[j]
    y2[j, 0] = pfY_array[j]
'''  
x2[:, 0] = x[:, Nt-1]
y2[:, 0] = y[:, Nt-1]

#T = np.linspace(0.01, 7, Nt)
T = 3.5

for n in range(0,Nt-1):
    energiaF = 0
    energiaP = 0
    Part = random.randint(0,N-1)
    varx = random.uniform(-D,D)
    vary = random.uniform(-D,D) 
    x2[:,n+1] = x2[:,n]
    y2[:,n+1] = y2[:,n]
    x2[Part, n+1] += varx
    y2[Part, n+1] += vary
    for k in range (N):
        if Part != k:
            energiaF += U(x2[Part,n+1], y2[Part, n+1], x2[k,n+1], y2[k,n+1])
            energiaP += U(x2[Part,n], y2[Part, n], x2[k,n], y2[k,n])
    DeltaE = (1/2)*(energiaF - energiaP)
    x2[Part,n+1] = x2[Part,n]
    y2[Part,n+1] = y2[Part,n]
    #print("ΔE=", DeltaE)
    if DeltaE > 0:
        u = random.uniform(0,1)
        #p = np.exp((-DeltaE)/(Kb*T[n]))
        p = np.exp((-DeltaE)/(Kb*T))
        if u > p:
            varx, vary = 0, 0
        
    x2[Part,n+1] += varx
    y2[Part,n+1] += vary


#Configurando a figura onde será feita a animação
dt = 0.01
L = 5
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(autoscale_on=False, xlim=(-L, L), ylim=(-L, L))
ax.set_aspect('equal')
    
# Objetos que receberão os dados (configurações) e texto a serem mostrados em cada quadro
confs, = ax.plot([], [], 'ko', ms=5, alpha=0.4)
texto = ax.text(0.05, 0.9, '', transform=ax.transAxes)

# Converte os dados recebidos em cada iteração em objetos a serem mostrados na figura
def animate(n):
    confs.set_data(x2[:,400*n], y2[:,400*n])
    texto.set_text('Temperatura = %.2f' % (T))
    return confs, texto

# Constrói a animação
ani = FuncAnimation(fig, animate, frames=Nt//400,interval=50)


plt.show()
ani.save("Agregados_atomicos_temperatura_3_5.gif")
