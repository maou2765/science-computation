import numpy as np
import matplotlib.pylab as plt

class const:
    def __init__(self,a,m,v,hbar):
        self.a=a
        self.m=m
        self.v=v
        self.hbar=hbar
        self.epi=2*m*a**2/hbar**2
        self.v*=self.epi
c=const(1.0e-11,9.11e-31,50*1.6022e-19,1.05457e-34)

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
    h=8/2400.
    for t0 in np.linspace(-4,4,2400):
            k1=h*Harmonic_virtual_4(sai,t0,E)
            k2=h*Harmonic_virtual_4(sai+k1/2.,t0+h/2.,E)
            k3=h*Harmonic_virtual_4(sai+k2/2.,t0+h/2.,E)
            k4=h*Harmonic_virtual_4(sai+k3,t0+h,E)
            sai+=(k1+2*k2+2*k3+k4)/6.

    return sai[0]

def wavefunction(E,boundary):
    good_phi=1.
    sai=np.array([0.,good_phi])
    h=np.abs(boundary[0])+np.abs(boundary[1])
    h/=6000.
    psi=[]
    for t0 in np.linspace(boundary[0],boundary[1],6000.):
        k1=h*Harmonic_virtual_2(sai,t0,E)
        k2=h*Harmonic_virtual_2(sai+k1/2.,t0+h/2.,E)
        k3=h*Harmonic_virtual_2(sai+k2/2.,t0+h/2.,E)
        k4=h*Harmonic_virtual_2(sai+k3,t0+h,E)
        sai+=(k1+2*k2+2*k3+k4)/6.
        psi.append(sai[0])
    tmp=0
    for psix in psi:
        tmp+=psix**2
    psi/=tmp**0.5
    return psi

def wavefunction_4(E,boundary):
    sai=np.array([0.,1.])
    h=np.abs(boundary[0])+np.abs(boundary[1])
    h/=2400.
    psi=[]
    for t0 in np.linspace(boundary[0],boundary[1],2400):
        k1=h*Harmonic_virtual_4(sai,t0,E)
        k2=h*Harmonic_virtual_4(sai+k1/2.,t0+h/2.,E)
        k3=h*Harmonic_virtual_4(sai+k2/2.,t0+h/2.,E)
        k4=h*Harmonic_virtual_4(sai+k3,t0+h,E)
        sai+=(k1+2*k2+2*k3+k4)/6.
        psi.append(sai[0])
        tmp=0
    for psix in psi:
        tmp+=psix**2
    psi/=tmp**0.5
    return psi

E_square=secant_method(f_2,3,[0.,1.5],1e-8)
print('-')
tmp_E_square=E_square.copy()
E_square=np.array(E_square)
E_square/=(c.epi*1.6e-19)
print(E_square)
print("Energy difference between the ground state and first excited state",E_square[1]-E_square[0])
print("Energy difference between the first excited state and secand excited state",E_square[2]-E_square[1])
print('-')
E_quad=[]
E_quad.append(secant_method(f_4,1,[-1.,0.2],1e-8)[0])
E_quad.append(secant_method(f_4,1,[1.,2.],1e-8)[0])
E_quad.append(secant_method(f_4,1,[0.8,1.],1e-8)[0])#E_quad.append(secant_method(f_4,1,[0.2,1.],1e-8)[0])
print(E_quad)
tmp_E_quad=E_quad.copy()
E_quad=np.array(E_quad)
E_quad/=(c.epi*1.6e-19)
print(E_quad)
print("Energy difference between the ground state and first excited state",E_quad[1]-E_quad[0])
print("Energy difference between the first excited state and secand excited state",E_quad[2]-E_quad[1])

t=np.linspace(-10,10,6000)
H_0=wavefunction(tmp_E_square[0]-0.01,[-10,10])
H_1=wavefunction(tmp_E_square[1],[-10,10])
H_2=wavefunction(tmp_E_square[2],[-10,10])
plt.plot(t,H_0,label="ground state energy-0.01")
plt.legend()
plt.ylabel('wavefunction')
plt.xlabel('x')
plt.title("not allowed energy level")
plt.show()
H_0=wavefunction(tmp_E_square[0],[-10,10])
plt.plot(t,H_0,label="ground state")
plt.legend()
plt.plot(t,H_1,label="First exited state")
plt.legend()
plt.plot(t,H_2,label="Second exited state")
plt.legend()
plt.ylabel('wavefunction')
plt.xlabel('x')
plt.title("wavefunction for V=v0x**2/a**2 (harmonic oscillator)")
plt.show()

H4_0=wavefunction_4(tmp_E_quad[0],[-4,4])
H4_1=wavefunction_4(tmp_E_quad[1],[-4,4])
H4_2=wavefunction_4(tmp_E_quad[2],[-4,4])
t=np.linspace(-4,4,2400)
plt.plot(t,H4_0,label="ground state")
plt.legend()
plt.plot(t,H4_1,label="First exited state")
plt.legend()
plt.plot(t,H4_2,label="Second exited state")
plt.legend()
plt.ylabel('wavefunction')
plt.xlabel('x')
plt.title("wavefunction for V=v0x**4/a**4(anharmonic oscillator)")
plt.show()

#symmetric calculation
def f_2_sym_odd(E):
    sai=np.array([0.,1.])
    h=10/3000.
    for t0 in np.linspace(0,10,3000.):
            k1=h*Harmonic_virtual_2(sai,t0,E)
            k2=h*Harmonic_virtual_2(sai+k1/2.,t0+h/2.,E)
            k3=h*Harmonic_virtual_2(sai+k2/2.,t0+h/2.,E)
            k4=h*Harmonic_virtual_2(sai+k3,t0+h,E)
            sai+=(k1+2*k2+2*k3+k4)/6.
    return sai[0]
def f_2_sym_even(E):
    sai=np.array([1.,0.])
    h=10/3000.
    for t0 in np.linspace(0,10,3000.):
            k1=h*Harmonic_virtual_2(sai,t0,E)
            k2=h*Harmonic_virtual_2(sai+k1/2.,t0+h/2.,E)
            k3=h*Harmonic_virtual_2(sai+k2/2.,t0+h/2.,E)
            k4=h*Harmonic_virtual_2(sai+k3,t0+h,E)
            sai+=(k1+2*k2+2*k3+k4)/6.
    return sai[0]
def wavefunction_sym_odd(E,boundary):
    good_phi=1.
    sai=np.array([0.,good_phi])
    h=np.abs(boundary[0])+np.abs(boundary[1])
    h/=3000.
    psi=[]
    for t0 in np.linspace(boundary[0],boundary[1],3000.):
        k1=h*Harmonic_virtual_2(sai,t0,E)
        k2=h*Harmonic_virtual_2(sai+k1/2.,t0+h/2.,E)
        k3=h*Harmonic_virtual_2(sai+k2/2.,t0+h/2.,E)
        k4=h*Harmonic_virtual_2(sai+k3,t0+h,E)
        sai+=(k1+2*k2+2*k3+k4)/6.
        psi.append(sai[0])
    tmp=0
    for psix in psi:
        tmp+=psix**2
    psi/=(2*tmp)**0.5
    return psi
def wavefunction_sym_even(E,boundary):
    sai=np.array([1.,0.])
    h=np.abs(boundary[0])+np.abs(boundary[1])
    h/=3000.
    psi=[]
    for t0 in np.linspace(boundary[0],boundary[1],3000.):
        k1=h*Harmonic_virtual_2(sai,t0,E)
        k2=h*Harmonic_virtual_2(sai+k1/2.,t0+h/2.,E)
        k3=h*Harmonic_virtual_2(sai+k2/2.,t0+h/2.,E)
        k4=h*Harmonic_virtual_2(sai+k3,t0+h,E)
        sai+=(k1+2*k2+2*k3+k4)/6.
        psi.append(sai[0])
    tmp=0
    for psix in psi:
        tmp+=psix**2
    psi/=(2*tmp)**0.5
    return psi
print('-')
E_square=[]
E_square_even_0=secant_method(f_2_sym_even,1,[0.0,0.4],1e-8)
E_square_even_1=secant_method(f_2_sym_even,1,[1.1,1.5],1e-8)
E_square_odd=secant_method(f_2_sym_odd,1,[0.9,1.1],1e-8)
E_square.append(E_square_even_0[0])
E_square.append(E_square_odd[0])
E_square.append(E_square_even_1[0])
print(E_square)
tmp_E_square=E_square.copy()
E_square=np.array(E_square)
E_square/=(c.epi*1.6e-19)
print(E_square)
print(1e-8/(c.epi*1.6e-19))
print("Energy difference between the ground state and first excited state",E_square[1]-E_square[0])
print("Energy difference between the first excited state and secand excited state",E_square[2]-E_square[1])
print('-')
t=np.linspace(0,10,3000)
H_0=wavefunction_sym_even(tmp_E_square[0],[0,10])
H_1=wavefunction_sym_odd(tmp_E_square[1],[0,10])
H_2=wavefunction_sym_even(tmp_E_square[2],[0,10])
plt.plot(t,H_0,label="ground state")
plt.legend()
plt.plot(t,H_1,label="First exited state")
plt.legend()
plt.plot(t,H_2,label="Second exited state")
plt.legend()
plt.ylabel('wavefunction')
plt.xlabel('x')
plt.title("wavefunction for V=v0x**2/a**2 (harmonic oscillator)")
plt.show()
