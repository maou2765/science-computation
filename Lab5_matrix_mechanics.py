import numpy as np
import matplotlib.pylab  as plt
import integration as ig
def Hmn(m,n):
    if (m!=n and (m%2)==(n%2)):
        return 0
    elif(m!=n and (m%2)!=(n%2)):
        return np.negative(80*(1.6022*10e-20)*m*n/(((m**2-n**2)**2)*(np.math.pi)**2))
    elif m==n:
        return 5*(1.6022e-19)+((n**2)*(np.math.pi**2)*(1.054571*10e-35)**2)/(50*9.1094*10e-52)
def H(order):
    H_list=[[Hmn(m,n) for n in range(1,order+1)] for m in range(1,order+1)]
    return np.array(np.matrix(H_list),subok=True)
E_10=np.linalg.eigvals(H(10))
E_100,v=np.linalg.eigh(H(100))
print(E_10[0]/(1.6022*10e-20))
print('10x10')
for i in E_10:
    print(i/(1.6022*10e-20))
i=0
print('-'*10)
print('100x100')
for j in E_100:
    print(j/(1.6022*10e-20))
    i+=1
    if i>=10:
        break
error=0
for i in range(0,10):
    error+=(E_10[i]-E_100[i])/2
error/=10
print(error)
print('-'*10)
print(np.linalg.norm(v[:,0]-v[0].T))
def wavefunction_n(v,n):
    return v[n].T
def prob_density(x,v,n):
    sai_n=0
    for i in range(0,100):
        sai_n+=wavefunction_n(v,n)[i]*np.math.sin(np.math.pi*(i+1)*x/5e-10)
    return sai_n
prob_sum=ig.numerical_method(ig.simpson,1,lambda x:prob_density(x,v,0)**2,lambda x,y: (x-y)/15.0,0,5e-10,100,10e-7)
prob_sum.adaptive_cal(0,5e-10,10e-8,0)
norm=prob_sum.latest_res()
v[0]/=np.sqrt(norm)
prob_sum=ig.numerical_method(ig.simpson,1,lambda x:prob_density(x,v,1)**2,lambda x,y: (x-y)/15.0,0,5e-10,100,10e-7)
prob_sum.adaptive_cal(0,5e-10,10e-8,0)
norm=prob_sum.latest_res()
v[1]/=np.sqrt(norm)
prob_sum=ig.numerical_method(ig.simpson,1,lambda x:prob_density(x,v,2)**2,lambda x,y: (x-y)/15.0,0,5e-10,100,10e-7)
prob_sum.adaptive_cal(0,5e-10,10e-8,0)
norm=prob_sum.latest_res()
v[2]/=np.sqrt(norm)
def prob_density_n(x,v,n):
    return np.asscalar(prob_density(x,v,n)**2)
x_axis=[x*5e-10/100 for x in range(0,101)]
y_1=[prob_density_n(x*5e-10/100,v,0) for x in range(0,101)]
y_2=[prob_density_n(x*5e-10/100,v,1) for x in range(0,101)]
y_3=[prob_density_n(x*5e-10/100,v,2) for x in range(0,101)]

plt.plot(x_axis,y_1,'b-',label=u'Ground state')
plt.legend()
plt.plot(x_axis,y_2,'r-',label=u'First excited state')
plt.legend()
plt.plot(x_axis,y_3,'g-',label=u'Second excited state')
plt.legend()

plt.ylabel("Probabilty density")
plt.xlabel('x')
plt.show()
