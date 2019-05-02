import numpy as np
import matplotlib.pylab as plt
def f_c(c):
    def f(x):
        return 1-np.exp(-1*c*x)
    return f
def relaxation(f,x,epi):
    i=1
    tmp=[x,f(x),f(f(x))]
    error=(tmp[1]-tmp[2])/(1-1/((tmp[2]-tmp[1])/(tmp[1]-tmp[0])))
    while np.abs(error)>epi:
        tmp[0]=tmp[1]
        tmp[1]=tmp[2]
        tmp[2]=f(tmp[2])
        i+=1
        error=(tmp[1]-tmp[2])/(1-1/((tmp[2]-tmp[1])/(tmp[1]-tmp[0])))
    return [tmp[2],i]
def over_relaxation(f,x,epi,omega):
    i=1
    tmp=[x,0,0]
    tmp[1]=f(x)*(1+omega)-omega*x
    tmp[2]=f(tmp[1])*(1+omega)-omega*tmp[1]
    error=(tmp[1]-tmp[2])/(1-1/((1+omega)*(tmp[2]-tmp[1])/(tmp[1]-tmp[0])-omega))
    while np.abs(error)>epi:
        tmp[0]=tmp[1]
        tmp[1]=tmp[2]
        tmp[2]=f(tmp[2])*(1+omega)-omega*tmp[2]
        i+=1
        error=(tmp[1]-tmp[2])/(1-1/((1+omega)*(tmp[2]-tmp[1])/(tmp[1]-tmp[0])-omega))
    return [tmp[2],i]
'''
root=relaxation(f_c(2),0.5,10e-7)
print(root[0],' with number of step:',root[1])

x_axis=np.linspace(0,3,300)
y=[relaxation(f_c(c),c/2.,10e-7)[0] for c in x_axis]
plt.plot(x_axis,y)
plt.show()
omega_list=[0,0.5,0.8,0.82,0.9,1,-0.1,-0.5,-0.8]
for omega in omega_list:
    root=over_relaxation(f_c(2),0.5,10e-7,omega)
    print(root,' omega=',omega)
'''
