import numpy as np
import matplotlib.pylab as plt

def runge_kutta_4(f,x0,t0,h):
    k1=h*f(x0,t0)
    k2=h*f(x0+k1/2.,t0+h/2.)
    k3=h*f(x0+k2/2.,t0+h/2.)
    k4=h*f(x0+k3,t0+h/2)
    return x0+(k1+2*k2+2*k3+k4)/6.
def runge_kutta_4_multivariable(f_list,r0,t0,h):
    for i in range(r0.shape[0]):
        k1=h*f_list[i](r0,t0)
        k1/=2.
        k2=h*f_list[i](r0+k1,t0+h/2.)
        k2/=2.
        k3=h*f_list[i](r0+k2,t0+h/2.)
        k4=h*f_list[i](r0+k3,t0+h)
        r0[i]+=(k1+2*k2+2*k3+k4)/6.
    return r0

def LPF(RC):
    def LPF_RC(x,t):
        def vin(t):
            if int((t//0.5)%2):
                return -1
            else:
                return 1
        return (vin(t)-x)/RC
    return LPF_RC
'''
f_list=[lambda r,t:10*(r[1]-r[0]),lambda r,t:28*r[0]-r[1]-r[0]*r[2],lambda r,t:r[0]*r[1]-8*r[2]/3.]

t_axis=np.linspace(0,10,10000)
y_1=[0]
for t in t_axis:
    y_1.append(runge_kutta_4(LPF(1),y_1[-1],t,0.001))
y_1.pop()
plt.plot(t_axis,y_1)
plt.show()
y_01=[0]
for t in t_axis:
    y_01.append(runge_kutta_4(LPF(0.1),y_01[-1],t,0.001))
y_01.pop()
plt.plot(t_axis,y_01)
plt.show()
y_001=[0]
for t in t_axis:
    y_001.append(runge_kutta_4(LPF(0.01),y_001[-1],t,0.001))
y_001.pop()
plt.plot(t_axis,y_001)
plt.show()
t_axis=np.linspace(0,50,5000)
r=[[0.,1.,0.]]
for t in t_axis:
    r.append(runge_kutta_4_multivariable(f_list,np.array(r[-1]),t,0.01))
y_axis=[]
for i in r:
    y_axis.append(i[1])
y_axis.pop()
plt.plot(t_axis,y_axis)
plt.show()
x_axis=[]
for i in r:
    x_axis.append(i[0])
x_axis.pop()
z_axis=[]
for i in r:
    z_axis.append(i[2])
z_axis.pop()
plt.plot(x_axis,z_axis)
plt.show()
'''
