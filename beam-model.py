# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 18:36:49 2019

@author: Dmitry

Данный файл позволяет выполнить расчёт теории взаимодействия проводящей плёнки с излучением в оптическом волноводе.
В качестве материалов рассматривается нитрид ниобия NbN или титан Ti (сейчас счёт для него, в NbN сложнее из-за мнимой части n,
смотри файл NbN-beam-model в соответствующей папке) в комбинации с волноводом в ниобате лития Ti:LiNbO3.
Изучена возможность создания плазмон поляритонных волн, также рассчитаны все условия для этого в случае как сплошного слоя, так и дифракционной решётки.
Всё это представимо в виде расчётов коэффициентов, однако переход от них к собственной плазмон-поляритонам описан в дипломной работе.
Рисунки переделаны для максимальной наглядности для использования в отчётах по выполнению задания.
Логика построения формул и обозначение величин объясняются в дипломной работе.
Убраны расчёты интеграла через scipy, можешь найти это в файле NIR.py, там показано, что в моём случае использовать их нецелесообразно.
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as s
import cmath as cm

#первый блок данных и констант
pi=np.pi
n1=3.6799
n2=4.5327
nT=n1+1j*n2
Rw=1
W1=7*10**-6
W2=10*10**-6
fi=pi/2-0.03
teta=pi/2-fi
n=2.142
la=1550*10**-9
c=3*10**8
w=2*pi*c/la
nv=1

#расчёт коэффициентов отражения для ТЕ и ТМ мод
r12E=(n*np.cos(fi)-n1*np.cos(teta)-1j*n2*np.cos(teta))/(n*np.cos(fi)+n1*np.cos(teta)+1j*n2*np.cos(teta))
r23E=(n1*np.cos(teta)+1j*n2*np.cos(teta)-(1-n1**2*np.sin(teta))**(1/2))/(n1*np.cos(teta)+1j*n2*np.cos(teta)+(1-n1**2*np.sin(teta))**(1/2))

r12M=(n1*np.cos(fi)-n*np.cos(teta)+1j*n2*np.cos(fi))/(n1*np.cos(fi)+n*np.cos(teta)+1j*n2*np.cos(fi))
r23M=(np.cos(teta)-n1*((1-n1**2*np.sin(teta))**(1/2))-1j*n2*((1-n1**2*np.sin(teta))**(1/2)))/(np.cos(teta)+n1*((1-n1**2*np.sin(teta))**(1/2))+1j*n2*((1-n1**2*np.sin(teta))**(1/2)))

d=np.arange(0,10**-6,10**-9) #рассматриваемый диапазон толщины слоя металла
delta=4*pi*d*nT*np.cos(teta)/la #набег фазы
RmE=(abs((r12E+r23E*np.exp(1j*delta))/(1+r12E*r23E*np.exp(1j*delta))))**2
RmM=(abs((r12M+r23M*np.exp(1j*delta))/(1+r12M*r23M*np.exp(1j*delta))))**2

#расчёт коэффициентов альфа (амплитудная часть), извлекание мнимой части
alphaE=(1-RmE)/(2*W1*np.tan(fi))
alphaM=(1-RmM)/(2*W1*np.tan(fi))
alphaE2=(1-RmE)/(2*W2*np.tan(fi))
alphaM2=(1-RmM)/(2*W2*np.tan(fi))
alEim=alphaE.imag
alMim=alphaM.imag
alEim2=alphaE2.imag
alMim2=alphaM2.imag

pog1=10*np.log10(np.exp(alphaE*10**-3))
pog2=10*np.log10(np.exp(alphaM*10**-3))

d1=d*10**6 #переназначение для удобства чтения графиков
print("_________________________________________")
print("Коэффициент r12 для TE-моды",round(r12E,2))
print("Коэффициент r23 для TE-моды",round(r23E,2))
print("Коэффициент r12 для TM-моды",round(r12M,2))
print("Коэффициент r23 для TM-моды",round(r23M,2))

fig, (ax1,ax2,ax3,ax4) = plt.subplots(nrows=4,ncols=1, figsize=(8,25))
ax1.plot(d1,RmE,d1,RmM)
ax1.axis([-10**-2, 3*10**-1, 0.1, 1.1])
ax1.set_xlabel("Толщина слоя металла, мкм")
ax1.set_ylabel("Коэффициент Rm")
ax1.grid()
ax1.set_title("Коэффициент Rm для TE (синий) и TM (оранжевый) мод")

ax2.plot(d1,alphaE,d1,alphaE2)
ax2.axis([-10**-2, 3*10**-1, 0, 700])
ax2.set_xlabel("Толщина слоя металла, мкм")
ax2.set_ylabel("Коэффициент alpha для TE-моды")
ax2.grid()
ax2.set_title("Графики alpha для TE, синий W=7мкм, оранжевый W=10мкм")

ax3.plot(d1,alphaM,d1,alphaM2)
ax3.axis([-10**-2, 3*10**-1, 0, 600])
ax3.set_xlabel("Толщина слоя металла, мкм")
ax3.set_ylabel("Коэффициент alpha для TM-моды")
ax3.grid()
ax3.set_title("Графики alpha для TM, синий W=7мкм, оранжевый W=10мкм")

ax4.plot(d1,pog1,d1,pog2)
ax4.axis([-10**-2, 1.5*10**-1, 0, 3])
ax4.set_xlabel("Толщина слоя металла, мкм")
ax4.set_ylabel("Погонные потери, дБ/мм")
ax4.grid()
ax4.axvline(x=0.008,ymin=0,ymax=3,color='red')
ax4.axvline(x=0.012,ymin=0,ymax=3,color='red')
ax4.set_title("Погонные потери для TE (синий) и TM (оранжевый) мод")
plt.show()
dmax=44*10**-9

pog10nmE=10*np.log10(np.exp(alphaE[10]*10**-3))
pog10nmM=10*np.log10(np.exp(alphaM[10]*10**-3))
print("_________________________________________________")
print("Погонные потери для TE на 10нм металла",round(pog10nmE,2),"дБ/мм")
print("Погонные потери для TM на 10нм металла",round(pog10nmM,2),"дБ/мм")

#второй блок данных и констант
e=1.602*10**-19
me=9.11*10**-31
N=12.475*10**29
eps0=8.85*10**-12
wp=((N*e**2)/(eps0*me))**(1/2)
kp=wp*c #плазменное волновое число
epsr=n1**2-n2**2
epsi=2*n1*n2
epsTi=epsr-1j*epsi
epsv=1
eps1=epsv
eps3=n**2
wsp1=wp*((1/2)**(1/2)) #Sp частота первого слоя
wsp3=wp*((1/(1+eps3))**(1/2)) #Sp частота второго слоя
wx=np.arange(1,6*10**15,10**12)
a=0.147*10**-9
T=300
k=1.38*10**-23

#определение дополнительных величин
wsp3=np.ones(wx.size)*wsp3 #прямая Sp1 частоты, для графиков
wsp1=np.ones(wx.size)*wsp1 #прямая Sp3 частоты, для графиков
war=np.ones(wx.size)*w #прямая частоты, для графиков
epsm=1-(wp/wx)**2
k0=wx/c
tau=(me**(1/2))/(N*pi*a**2*(3*k*T)**(1/2)) #время релаксации
kx=wx/c
kx2=wx/(c/n)
wpp=(wp**2+c**2*kx**2)**(1/2)
lax=2*pi*c/wx
lapp=2*pi*c/wpp
lasp3=2*pi*c/wsp3
lasp1=2*pi*c/wsp1
laar=2*pi*c/war

#пересчёт волновых чисел через эпсилоны, для графиков
kxn1=((wx/c)**2/(epsm*eps1/(epsm+eps1)))**(1/2)
kxn3=((wx/c)**2/(epsm*eps3/(epsm+eps3)))**(1/2)
w3=c*kx*((epsm+eps3)/(epsm*eps3))**(1/2)
w1=c*kx*((epsm+eps1)/(epsm*eps1))**(1/2)
la1=2*pi*c/w1
la3=2*pi*c/w3

fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(7,5))
ax1.plot(kx,wx,kx2,wx)
#plt.plot(kx,wpp)
#plt.plot(kx,wsp3,'--',kx,wsp1,'--')
ax1.plot(kx,w3,kx,w1,kx,war)
ax1.grid()
ax1.axis([-10**6,1.5*10**7,-5*10**14,2*10**15])
#plt.axis([-10**4,3*10**5,-5*10**12,10*10**13])
ax1.set_ylabel("Частота, w, рад/с")
ax1.set_xlabel("Фазовый коэффициент, kx, 1/м")
plt.show()

#ищу точки пересечения дисперсионных кривых с прямой исследуемой длины волны
r1 = np.argwhere(np.diff(np.sign(war - w1))).flatten()
r3 = np.argwhere(np.diff(np.sign(war - w3))).flatten()
print("Точки пересечения по kx:",round(kx[r1][0]/10**6,1),"* 10^6;",round(kx[r3][0]/10**6,1),"* 10^6 1/м")
#print("Точка пересечения по w:",round(w1[r1][0]/10**15,3),"Трад/с")
print("Точки пересечения по la:",round(la1[r1][0]*10**6,2),"мкм")

dkx=abs(kx[r1][0]-kx[r3][0]) #разница в волновом числе
lap=2*pi/dkx #пересчёт разницы в длину волны
b=lap/2
dim=10**-3
Nr=dim/lap
hi=np.arange(0,2*pi,pi/10000)
I0=1*(((np.sin((np.sin(fi)-np.sin(hi))*pi*b/la)))/((np.sin(fi)-np.sin(hi))*pi*b/la))**2*(((np.sin((np.sin(fi)-np.sin(hi))*pi*Nr*lap/la)))/((np.sin(fi)-np.sin(hi))*pi*lap/la))**2
grN=-0.05*np.max(I0)
grV=1.1*np.max(I0)

#def U(psi): return (((np.sin((np.sin(fi)-np.sin(psi))*pi*b/la)))/((np.sin(fi)-np.sin(psi))*pi*b/la))**2*(((np.sin((np.sin(fi)-np.sin(psi))*pi*Nr*lap/la)))/((np.sin(fi)-np.sin(psi))*pi*lap/la))**2

print("_______________________________________________")
print("Период решётки",round(lap*10**6,3),"мкм")
nr=1
Q=2*pi*la*dmax/(nr*lap*np.cos(fi)) #эффективная толщина решётки

fi2=cm.asin(n*np.sin(fi)/nv) #расчёт угла распространения луча в основной среде
RperpTE=(((cm.sin(fi2-fi))**2)/((cm.sin(fi2+fi))**2)) #коэффициент отражения для параллельной поляризации (TE)
RparalTM=(((cm.tan(fi2-fi))**2)/((cm.tan(fi2+fi))**2)) #коэффициент отражения для перпендикулярной поляризации (TM)

#расчёт фазовых сдвигов для 10нм толщины материала (металла или сверхпроводника), формулы согласно дипломке
RMM=RmM[10]**(1/2)
RME=RmE[10]**(1/2)
RTM=abs(RparalTM)
RTE=abs(RperpTE)
eiE=RperpTE/RTE
eiM=RparalTM/RTM
E1=np.arctan(eiE.imag/eiE.real)
M1=np.arctan(eiM.imag/eiM.real)
print("Отражается TE",round(RTE*100,1),"%")
print("Отражается TM",round(RTM*100,1),"%")
YM=((r12M+r23M*np.exp(1j*delta))/(1+r12M*r23M*np.exp(1j*delta)))**2
YE=((r12E+r23E*np.exp(1j*delta))/(1+r12E*r23E*np.exp(1j*delta)))**2
AM=abs(YM)
AE=abs(YE)
eiFM=YM/AM
eiFE=YE/AE
FE=np.arctan(eiFE.imag/eiFE.real)
FM=np.arctan(eiFM.imag/eiFM.real)
dFM=FM-M1
dFE=FE-E1

#получаемый фазовый сдвиг в продольной компоненте волнового вектора на одиночной длине переотражения
dbetaM1=dFM/(2*W1*np.tan(fi))
dbetaM2=dFM/(2*W2*np.tan(fi))
dbetaE1=dFE/(2*W1*np.tan(fi))
dbetaE2=dFE/(2*W2*np.tan(fi))
 
#блок графиков, внимательно подбирай коэффициенты для графиков, лучше вообще отрубай их при смене входных данных
fig, (ax1, ax2) = plt.subplots(nrows=1,ncols=2,figsize=(17,5))
ax1.plot(d1,np.rad2deg(FE),d1,np.rad2deg(FM))
#ax1.plot(d1,np.rad2deg(dFE),d1,np.rad2deg(dFM))
ax1.grid()
ax1.axis([-0.01,0.3,-17,12])
ax1.set_ylabel("Фаза мнимой части, град")
ax1.set_xlabel("Толщина слоя металла, мкм")
ax1.set_title("Синий - TE, оранжевый - TM, зелёный - delta TE, красный - delta TM")
ax2.plot(d1,dbetaE1,d1,dbetaM1)
ax2.grid()
ax2.axis([-0.002,0.15,-600,400])
ax2.set_ylabel("Delta beta, 1/м")
ax2.set_xlabel("Толщина слоя металла, мкм")
ax2.axvline(x=0.008,ymin=0,ymax=3,color='red')
ax2.axvline(x=0.012,ymin=0,ymax=3,color='red')
#ax2.set_title("Синий и оранжевый - TM W=7 и 10мкм, зелёный и красный - TE W=7 и 10мкм")
plt.show()

print("Максимум fi",round(np.min(FM),3),"рад|",round(np.rad2deg(np.min(FM)),1),"град при толщине слоя металла",round(d1[np.argwhere(FM==np.min(FM))[0][0]],3),"мкм")
print("Значение на хвосте",round(FM[400],3),"рад|",round(np.rad2deg(FM[400]),1),"град")

#переводы для удобства представления на графике
dnM1=dbetaM1*la/(2*pi)*10**5
dnM2=dbetaM2*la/(2*pi)*10**5
dnE1=dbetaE1*la/(2*pi)*10**5
dnE2=dbetaE2*la/(2*pi)*10**5
fig, (ax2) = plt.subplots(nrows=1,ncols=1,figsize=(8,5))
ax2.plot(d1,dnM1,d1,dnM2,d1,dnE1,d1,dnE2)
ax2.grid()
ax2.axis([-0.01,0.3,-15,10])
ax2.set_ylabel("delta n, 10^-5")
ax2.set_xlabel("Толщина слоя металла, мкм")
ax2.set_title("Синий и оранжевый - TM W=7 и 10мкм, зелёный и красный - TE W=7 и 10мкм")
plt.show()

#%%
#блок расчёта дифракционной решётки, с целью создания плазмон-поляритонных волн не только в TM области

#рисование повторяющейся дельта-функции как вид функции пропускания
def S(x):
    for i in range(int(Nr)):
        if lap*i<=x and x<lap*i+b:
             return 1
        elif lap*i+b<=x and x<lap*i+2*b:
             return RMM**2

Nsh=30000 #число точек, выбрано довольно рандомно, главное, чтоб дельта функция выглядела нормально
sh=(dim-lap)/Nsh
full=np.linspace(0,dim-lap,Nsh,endpoint=False)
Sd=np.array([S(x) for x in full]) #перевод функции в точки
fSd=(np.fft.fft(Sd)) #Фурье-образ функции пропускания

nch=1 #стартовая точка, 0 лучше не ставить, так как там ошибка счёта
Tdim=lap
ffull=np.linspace(0,1/(2*Tdim),Nsh//2,endpoint=False)

full3=full*1000
fig, ((ax1),(ax2),(ax3)) = plt.subplots(nrows=3,ncols=1,figsize=(7,15))
ax1.plot(full3,Sd)
ax1.axis([0,2*10**-2,0.82,1.05])
ax1.set_title('Функция пропускания')
ax1.set_xlabel('$x*10^-3$ м')
ax1.set_ylabel('$T(x)$')
ax2.plot(ffull[nch:],fSd[nch:Nsh//2])
#ax2.set_ylim((-60,40))
ax2.set_title('Фурье-образ')
ax2.set_xlabel('$kx$ 1/м')
ax2.set_ylabel('$f(T(x))$')
ax3.plot(ffull[nch:],abs(fSd[nch:Nsh//2]))
#ax3.set_ylim((-400,250))
ax3.set_title('Модуль Фурье-образа')
ax3.set_xlabel('$kx$ 1/м')
ax3.set_ylabel('$|f(T(x))|$')
#ax4.plot(ffull[nch:],abs(fSd[nch:Nsh//2])**2)
#ax4.set_ylim((-1000,5000))
#ax4.set_title('Квадрат модуля Фурье-образа')
#ax4.set_xlabel('$kx$ 1/м')
#ax4.set_ylabel('$|f(T(x))|^2$')               не включай, если не нужен квадрат модуля Ф-обр, он не имеет смысла
plt.show()
print("Число точек при счёте =", Nsh)

aSd=abs(fSd)

I1=np.trapz((abs(fSd[1:int(Nr)*2]))) #расчёт интеграла через метод трапеции
I2=np.trapz((abs(fSd[1:Nsh//2])))
check=I1/I2*100
print("Отражается обратно на решётке",round(check,1),"%")
print("Проходит (диффракционная эффективность)",round(100-check,1),"%")
print("________________________________________________")
print("Отражается обратно при титане",round(RmM[10]*100,1),"%")
print("Проходит при титане",round(100-RmM[10]*100,1),"%")

#plt.plot(d1,F1,d1,pi-F2)
#plt.axis([-10**-2, 3*10**-1, -0.2, 3.3])
#plt.xlabel("Толщина слоя металла, мкм")
#plt.ylabel("Фазовый коэффициент")
#plt.grid()
#plt.show()
#print("Графики fi с разницей в pi")
#print("delta_fi",round(np.rad2deg(np.max(F2)-np.min(F2)),1),"град")












































#wx=np.arange(1,4*wsp,wsp/10)
#epsm=1-(wp/wx)**2
#kx=(wx*(epsm*epsv/(epsm+epsv))**(1/2))/c
#q=kx/kp
#wgr=wp*(1+1/(2*q**2)+(1+1/(4*q**4))**(1/2))**(-1/2)
#kxr=(wx*(epsr*epsv/(epsr+epsv))**(1/2))/c
#kxi=((wx*(epsr*epsv/(epsr+epsv))**(3/2))/c)*(epsi/(2*epsr**2))

#beta1r=(eps1)**(1/2)*k0*((epsr**2-epsr*eps1+eps1**2)/((abs(epsTi+eps1))**2))**(1/2)
#beta1i=(eps1**2*epsi*k0**2)/(beta1r*(abs(epsTi+eps1))**2*2)
#beta3r=(eps3)**(1/2)*k0*((epsr**2-epsr*eps3+eps3**2)/((abs(epsTi+eps3))**2))**(1/2)
#beta3i=(eps3**2*epsi*k0**2)/(beta1r*(abs(epsTi+eps3))**2*2)
#beta1=(beta1r**2+beta1i**2)**(1/2)
#beta3=(beta3r**2+beta3i**2)**(1/2)
#S11=(beta1**2-eps1*k0**2)**(1/2)
#S12=(beta1**2-epsTi*k0**2)**(1/2)
#S13=(beta1**2-eps3*k0**2)**(1/2)

#plt.plot(beta1r,wx,beta3r,wx)
#plt.axis([-0.2*10**10, 6*10**10, -0.1*10**-8, 10**-8])
#plt.xlabel("Толщина слоя металла, мкм")
#plt.ylabel("Погонные потери, дБ/мм")
#plt.grid()
#plt.show()
#plt.plot(beta1i,wx,beta3i,wx)
#plt.axis([-0.2*10**10, 6*10**10, -0.1*10**-8, 10**-8])
#plt.xlabel("Толщина слоя металла, мкм")
#plt.ylabel("Погонные потери, дБ/мм")
#plt.grid()
#plt.show()
#print("Погонные потери для TE (синий) и TM (оранжевый) мод")








