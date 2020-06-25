# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt 



# Shell and Tube heat exchanger

L = 60 # m, length of pipe
r1 = 0.1 # m, inner pipe diameter
r2 = 0.15 # m, outer pipe diameter
n = 100 # number of increments/nodes

m1 = 3 # kg/s, mass flow rate of fluid 1
Cp1 = 4180 # J/Kg*K, specific heat capacity of fluid 1 (water)
rho1 = 1250 # kg/m(^3), density of fluid 1 (crude oil) 

m2 = 5 # kg/s, mass flow rate of fluid 2
Cp2 = 4180 # J/kg*K, specific heat capacity of fluid 2 (water)
rho2 = 1000 # kg/m(^3), density of fluid 2 (water)

pi = 3.141592653589793
t_w = 2 # wall thickness,
Ac1 = pi*r1**2 # m^2, cross-sectional area of inner pipe
Ac2 = pi*(r2**2-r1**2) # cross sectional area of outer pipe

T1i = 400 # K, inlet temperature of pipe 1
T2i = 800 # K, surface temperature of pipe 2
T0 = 300 # K, initial temperature of fluid throughout the pipe

 
k = 18 # W/m*K, thermal conductivity of pipe shell side
Rflow = R2-dth # Flow radius of fluid 
dth = # fouling layer thickness 
hi = (k*Nu)/(2*Rflow) # tube side heat transfer co-efficient
Nu=0.0023*(Re**0.8)*(Pr**0.4) # Dittus Boelter Correlation for Nusselt number
Re = (rho1*v*d)/mu) # Reynolds number
Pr = 6 # Prandtl Number (arbitrary value)
U=(1/hi)+(Rfn)+(t_w/k)+(1/ho) # overall heat transfer co-efficient


dx= L/n # cell width

t_final = 1000 #s, simulation time 
dt= 1 # s, time step

x = np.linspace(dx/2, L-dx/2, n)

T1 = np.ones(n)*T0
T2 = np.ones(n)*T0

dT1dt = np.zeros(n)
dT2dt = np.zeros(n)

t = np.arange(0, t_final, dt)

for j in range(1,len(t)): 
    
    plt.clf()
    
    dT1dt[1:n] = (m1*Cp1*(T1[0:n-1]-T1[1:n])+U*2*pi*r1*dx*(T2[1:n]-T1[1:n]))/(rho1*Cp1*Ac1*dx) # energy balance from node 1 to 100 for fluid 1
    dT1dt[0] = (m1*Cp1*(T1i-T1[0])+U*2*pi*r1*dx*(T2[0]-T1[0]))/(rho1*Cp1*Ac1*dx) # boundary condition for fluid 1
    
    dT2dt[1:n] = (m2*Cp2*(T2[0:n-1]-T2[1:n])-U*2*pi*r1*dx*(T2[1:n]-T1[1:n]))/(rho2*Cp2*Ac2*dx) # energy balance from node 1 to 100 for fluid 2
    dT2dt[0] = (m2*Cp2*(T2i-T2[0])-U*2*pi*r1*dx*(T2[0]-T1[0]))/(rho2*Cp2*Ac2*dx) # boundary condition for fluid 2
    
    T1 = T1 + dT1dt*dt
    T2 = T2 + dT2dt*dt
    
    plt.figure(1)
    plt.plot(x, T1, color = 'blue', label='inside')
    plt.plot(x, T2, color = 'red', label='outside')
    plt.axis ([0, L, 298, 820])
    plt.xlabel('Distance (n)')
    plt.ylabel('Temperature (K)')
    plt.legend(loc = 'upper right')
    plt.show()
    plt.pause(0.005)
              
    
    
   
    
    
