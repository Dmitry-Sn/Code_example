# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 12:09:09 2018

@author: Dmitry
Этот код позволяет выполнить симуляцию распределения амплитуды поля по слоистой структуре, а также подсчитать процент ухода поля в слой сверхпроводника.
Доработка файла NIR.py на случай слоя кремния.
"""

import math as m
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as s

#первый блок данных
n11=5.23
n12=5.82
n1=n11+1j*n12
nv=1
nsi=3.5
#print(n1)
#n1=(n1+36*nsi)/37
#print(n1)
#n1=(n1+nv)/2
#print(n1)
n1=(n1+nsi)/2
n11s=(n11+nsi)/2
n2=2.216
n3=2.214
n0=1
la=1550*10**-9
d=5*10**-6
a=10*10**-9
pi=3.1415926
k=(2*pi)/la
alpha=1.50098
g=k*n0*m.sin(alpha)
A=5
B=1
l=d/2
d1=d/2+a

#расчёт поперечных волновых векторов
h=k*n2*m.sin(alpha)
v0=(h**2*k**2)**(1/2)
v3=(h**2-(k*n3)**2)**(1/2)
u1=((k*n1)**2-h**2)**(1/2)
u2=((k*n2)**2-h**2)**(1/2)

#разделил части для расчёта критических размеров слоёв, чтоб не нагромождать код в длину
fst=n2**2-n0**2
snd=nsi**2-n2**2
trd=(fst/snd)**(1/2)
fth=m.atan(trd)
ffth=2*pi*snd**(1/2)
sth=fth/ffth

#расчёт критических размеров слоёв (NbN, Si, Si+NbN соответственно)
hkr=la*m.atan((((n2**2-n0**2)/(n11**2-n2**2))**(1/2)))/(2*pi*((n11**2-n2**2)**(1/2)))
hkr2=la*sth
hkr3=la*m.atan((((n2**2-n0**2)/(n11s**2-n2**2))**(1/2)))/(2*pi*((n11s**2-n2**2)**(1/2)))

#расчёт точек для графика распределения амплитуды поля по заданной слоистой структуре
x=np.arange(-l+a/1000,0.000004,a/10)
x1=np.arange(l-a,d1+a,a/1000)
def Et(x):
    if d1 < x:
        return A*np.exp(-v0*x)
    elif (l <= x) and (x <= d1):
        return A*np.exp(1j*u1*x)+B*np.exp(-1j*u1*x)
    elif (-l < x) and (x < l):
        return A*np.exp(1j*u2*x)+B*np.exp(-1j*u2*x)
    elif (x <= -l):
        return A*np.exp(v3*x)
    else:
        return 0
   
#интеграл методом трапеции (грубая оценка)
def trapezoidal(f, a, b, n):
    h = (b-a)/float(n)
    I = f(a) + f(b)
    for k in range(1, n, 1):
        x = a + k*h
        I += 2*f(x)
    I *= h/2
    return I
    
y = np.vectorize(Et, otypes=[float])
fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(10,5))
ax1.axvline(x=l,ymin=0,ymax=3,color='red')
ax1.set_ylabel("Амплитуда поля, В/м")
ax1.set_xlabel("Координата, м")
ax1.axis([-0.000003,0.000004,-1, 6.5])
ax1.grid()
ax1.set_title("подложка                              сердцевина                                               плёнка              воздух      ")
#ax1.set_title("") 
#plt.axvline(x=-l,ymin=0,ymax=3,color='red')
ax1.axvline(x=d1,ymin=0,ymax=3,color='red')
ax1.axhline(y=0,xmin=-10,xmax=10,color='red')
ax1.plot(x,y(-x))
plt.show()

#расчёт интегралов
I1=trapezoidal(y,-l,d1,100)
I2=trapezoidal(y,l,d1,100)
I3=s.quad(y,-l,d1)
I4=s.quad(y,l,d1)
It1=I3[0]
It2=I4[0]
I5=0.3*a
Itog=I2*100/I1
Itog2=It2*100/It1
print("Усредненный показатель преломления",n1)
print("Процент в NbN (ч/з trapezoidal)",round(Itog,1),"%")
print("Процент в NbN (ч/з scipy.integrate.quad)",round(Itog2,1),"%")
print("Итоговая эффективность",round(Itog*0.3,1),"%")
print("Минимальный размер слоя NbN",round(hkr*10**9,1),"нм")
print("Минимальный размер слоя Si",round(hkr2*10**9,1),"нм")
print("Минимальный размер слоя Si+NbN",round(hkr3*10**9,1),"нм")
print("Возможно ли вообще распространение волны с длиной 1550нм в слое Si в 10нм?")
print(" ")
#не знаю, почему так выход с scipy интегралом, но это зависит видимот от шага, т.к. при увеличении слоя он становится нормальным