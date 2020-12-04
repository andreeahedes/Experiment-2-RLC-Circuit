# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 09:33:47 2020

@author: ah920
"""
import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt 

def function(omega,N, gamma,omega0):
     a=(N)/((omega**2*gamma**2)+(omega**2-omega0**2)**2)**(1/2)
     return a   
def Qfactor(fit):
    return max(function(frequency,*fit[0]))/np.std(function(frequency,*fit[0]), ddof=1)

L=1e-3
C=100e-9
omega0=15.91e3*(2*np.pi)
r1=1
r2=2
r3=3
    
r1f,r1V0,r1Vc,r1phase=np.loadtxt('Week 8/Resistance 1.csv', skiprows=1, delimiter=',', unpack=True)
r2f,r2V0,r2Vc,r2phase=np.loadtxt('Week 8/Resistance 2.csv', skiprows=1, delimiter=',', unpack=True)
r3f,r3V0,r3Vc,r3phase=np.loadtxt('Week 8/Resistance 3.csv', skiprows=1, delimiter=',', unpack=True)

r2fbun=np.delete(r2f,8)
r2Vcbun=np.delete(r2Vc,8)


r1Vc=r1Vc/(2*np.sqrt(2))
r2Vc=r2Vc/(2*np.sqrt(2))
r2Vcbun=r2Vcbun/(2*np.sqrt(2))
r3Vc=r3Vc/(2*np.sqrt(2))



r1f=2*np.pi*r1f
frequency=np.linspace(2*np.pi*6e3,2*np.pi*22e3,100000)
initialGuess=[10,1e4,omega0]
fit1=spo.curve_fit(function,r1f,r1Vc,p0=initialGuess, maxfev=10**6)
plt.scatter(r1f,r1Vc, label='Data Measured')
plt.plot(frequency, function(frequency,*fit1[0]),label='Fitting Curve')
plt.grid()
plt.legend()
plt.xlabel('Angular Frequency (s^-1)')
plt.ylabel('Voltage Across the Capacitor (V)')
plt.title('Resistance 1 Ohm')
plt.show()

r2fbun=2*np.pi*r2fbun
r2f=2*np.pi*r2f
frequency=np.linspace(2*np.pi*6e3,2*np.pi*23e3,1000000)
initialGuess=[10,1e4,omega0]
fit2=spo.curve_fit(function,r2fbun,r2Vcbun,p0=initialGuess, maxfev=10**6)
plt.scatter(r2fbun,r2Vcbun,label='Data measured')
plt.errorbar(r2f[8],r2Vc[8],yerr=3*r2Vc[8]/100, xerr=0.01e4, capsize=2, color='orange', label='Erroneous measurment')
plt.plot(frequency, function(frequency,*fit2[0]),label='Fitting curve')
plt.grid()
plt.legend()
plt.xlabel('Angular Frequency (s^-1)')
plt.ylabel('Voltage Across the Capacitor (V)')
plt.title('Resistance 2 Ohm')
plt.show()

r3f=2*np.pi*r3f
frequency=np.linspace(2*np.pi*6e3,2*np.pi*22e3,100000)
initialGuess=[10,1e4,omega0]
fit3=spo.curve_fit(function,r3f,r3Vc,p0=initialGuess, maxfev=10**6)
plt.scatter(r3f,r3Vc, label='Data Measured')
plt.plot(frequency, function(frequency,*fit3[0]), label='Fitting Curve')
plt.grid()
plt.legend()
plt.xlabel('Angular Frequency (s^-1)')
plt.ylabel('Voltage Across the Capacitor (V)')
plt.title('Resistance 3 Ohm')
plt.show()

print('The quality factors are Q1=',Qfactor(fit1),'Q2=',  Qfactor(fit2),'and Q3=', Qfactor(fit3))








