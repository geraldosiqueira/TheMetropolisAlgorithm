"""
Created on Fri Sep  9 10:33:31 2022

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


T = 0
Nt = 600000
N = 37
D = 0.002*a
k = 1

L = 2

x = np.zeros( [N, Nt] )
y = np.zeros( [N, Nt] )
x[0,0] = random.uniform(-L,L)
y[0,0] = random.uniform(-L,L)    
for l in range(1,N-1):
    for g in range(N):
        if l != g:
            while ((x[l,0]-x[g,0])*(x[l,0]-x[g,0]) + (y[l,0]-y[g,0])*(y[l,0]-y[g,0]))**2 < a :
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
    
for k in range(0,N-1):
    f = open("posicaofinalX.txt", "a")
    g = open("posicaofinalY.txt", 'a')
    f.write(str(x[k, Nt-1]) + ";")   
    g.write(str(y[k, Nt-1]) + ";")
    f.close()
    g.close()
    
dt = 0.01

#Configurando a figura onde será feita a animação
L = 5
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(autoscale_on=False, xlim=(-L, L), ylim=(-L, L))
ax.set_aspect('equal')
    
# Objetos que receberão os dados (configurações) e texto a serem mostrados em cada quadro
confs, = ax.plot([], [], 'ko', ms=5, alpha=0.4)
texto = ax.text(0.05, 0.9, '', transform=ax.transAxes)

# Converte os dados recebidos em cada iteração em objetos a serem mostrados na figura
def animate(n):
    confs.set_data(x[:,100*n], y[:,200*n])
    texto.set_text('Tempo = %.2f' % (n*dt))
    return confs, texto

# Constrói a animação
ani = FuncAnimation(fig, animate, frames=Nt//200,interval=50)


plt.show()
#ani.save("animacao_agregados_atomicos_37.gif")
