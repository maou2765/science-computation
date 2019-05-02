import numpy as np
import matplotlib.pylab as plt
import Lab8_rungekutta as runku

class const:
    def __init__(self,a,m,v,hbar):
        self.a=a
        self.m=m
        self.v=v
        self.hbar=hbar
        self.epi=2*m*a**2/hbar**2
        self.v*=self.epi
c=const(1.0e-11,9.11e-31,50*1.6022e-19,1.05457e-34)
def Harmonic_virtual_2(r,t,E):
    return np.array([r[1],(c.v*(t**2.)-E)*r[0]])

def f_2(E):
    sai=np.array([0.,1.])
    h=20/6000.
    for t0 in np.linspace(-10,10,6000.):
            k1=h*Harmonic_virtual_2(sai,t0,E)
            k2=h*Harmonic_virtual_2(sai+k1/2.,t0+h/2.,E)
            k3=h*Harmonic_virtual_2(sai+k2/2.,t0+h/2.,E)
            k4=h*Harmonic_virtual_2(sai+k3,t0+h,E)
            sai+=(k1+2*k2+2*k3+k4)/6.
    return sai[0]

def Harmonic_virtual_4(r,t,E):
    return np.array([r[1],(c.v*(t**4.)-E)*r[0]])

def f_4(E):
    sai=np.array([0.,1.])
    h=16/5000.
    for t0 in np.linspace(-8,8,5000.):
            k1=h*Harmonic_virtual_4(sai,t0,E)
            k2=h*Harmonic_virtual_4(sai+k1/2.,t0+h/2.,E)
            k3=h*Harmonic_virtual_4(sai+k2/2.,t0+h/2.,E)
            k4=h*Harmonic_virtual_4(sai+k3,t0+h,E)
            sai+=(k1+2*k2+2*k3+k4)/6.
    return sai[0]
def MD(order):
    mat=' -2. 1'+(order-2)*' 0'+' ;'
    mat+=' 1 -2 1 '+(order-3)*' 0'+' ;'
    i=3
    while i<order:
        mat+=(i-2)*' 0 '+'1 -2 1'+(order-i-1)*' 0'+' ;'
        i+=1
    mat+=' 0 '*(order-2)+'1 -2'
    mat=np.matrix(mat)
    return mat

def matrix_relaxation(f,epi,x0):
    x0,tmp=np.matmul(f,x0),x0
    err=np.linalg.norm(x0-tmp)
    while err>epi:
        x0,tmp=np.matmul(f,x0),x0
        err=np.linalg.norm(x0-tmp)
    return x0

def wavefunction(E,boundary):
    D=MD(1000)
    t_axis=np.linspace(boundary[0],boundary[1],1000)
    i=0
    V=np.zeros([1000,1000])
    for t in t_axis:
        V[i,i]=1./(E-c.v*t**2)
        i+=1
    D=np.matmul(V,D)
    sai=matrix_relaxation(D,1e-8,np.matrix(1000*' 3.').T)
    norm=0
    i=0
    while i<1000:
        norm+=sai[i,0]**2
        i+=1
    sai/=norm**0.5
    return sai

def wavefunction_4(E,boundary):
    D=MD(1000)
    t_axis=np.linspace(boundary[0],boundary[1],1000)
    i=0
    for t in t_axis:
        D[i]/=(E-c.v*t**4)
        i+=1
    sai=matrix_relaxation(D,1e-8,np.matrix(1000*' 2.').T)
    norm=0
    for i in sai:
        norm+=np.linalg.norm(i)**2
    sai/=norm**0.5
    psi=[0]
    for i in sai:
        psi.append(i)
    psi.append(0)
    return sai

def secant_method(f,total_no_of_root,interval,epi):
    offset=(interval[1]-interval[0])/total_no_of_root
    E1,E2=interval[0],interval[0]+offset
    tmp=[0,f(E1)]
    no_of_root=0
    roots=[]
    while no_of_root<total_no_of_root:
        while np.abs(E2-E1)>epi:
            tmp[0],tmp[1]=tmp[1],f(E2)
            df=E2-E1
            df/=(tmp[1]-tmp[0])
            E1,E2=E2,E2-tmp[1]*df
        no_of_root+=1
        roots.append(E2)
        interval[0]+=offset
        E1,E2=interval[0],interval[0]+offset
        tmp=[roots[-1],f(E1)]
    return roots

E_square=secant_method(f_2,3,[0.,1.5],1e-8)
'''
arr=[]
for ans in E:
    if ans!="NA":
        arr.append(ans)
'''
print('-')
print(E_square)
print(f_2(E_square[0]))
tmp_E_square=E_square.copy()
E_square=np.array(E_square)
E_square/=(c.epi*1.6e-19)
print(E_square)
print(1e-8/(c.epi*1.6e-19))
print("Energy difference between the ground state and first excited state",E_square[1]-E_square[0])
print("Energy difference between the first excited state and secand excited state",E_square[2]-E_square[1])
print('-')
E_quad=secant_method(f_4,1,[0.3,0.5],1e-8)
E_quad.append(secant_method(f_4,1,[0.8,1.],1e-8)[0])
E_quad.append(secant_method(f_4,1,[9.,10.],1e-8)[0])
print(E_quad)
tmp_E_quad=E_quad.copy()
E_quad=np.array(E_quad)
E_quad/=(c.epi*1.6e-19)
print(E_quad)
t=np.linspace(-5,5,1000)
H_0=wavefunction_4(tmp_E_quad[0],[-4,4])
H_1=wavefunction_4(tmp_E_quad[1],[-4,4])
H_2=wavefunction_4(tmp_E_quad[2],[-4,4])
plt.plot(t,H_0)
plt.plot(t,H_2)
plt.plot(t,H_1)
plt.show()
t=np.linspace(-10,10,1000)
H2_0=wavefunction(tmp_E_square[0],[-10,10])
H2_1=wavefunction(tmp_E_square[1],[-10,10])
H2_2=wavefunction(tmp_E_square[2],[-10,10])
plt.plot(t,H2_0)
plt.plot(t,H2_1)
plt.plot(t,H2_2)
plt.show()
