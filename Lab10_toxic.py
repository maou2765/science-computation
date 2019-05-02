import numpy as np
import matplotlib.pylab as plt
phi=np.zeros([101,101],dtype=np.float32)
roll=np.zeros([101,101])
roll[20:41,20:41]=-1./1e4
roll[60:81,60:81]=1./1e4
epi=1e-6
dphi=phi.copy()
delta=epi+1.
while delta>epi:
    dphi[1:100,1:100]=(phi[2:101,1:100]+phi[0:99,1:100]+phi[1:100,2:101]+phi[1:100,0:99]+roll[1:100,1:100])/4.
    delta=np.max(np.abs(phi-dphi))
    print(delta)
    phi,dphi=dphi,phi

x=np.linspace(0,1,100)
y=np.linspace(0,1,100)
X,Y=np.meshgrid(x,y)
dy,dx=np.gradient(-phi)
skip=(slice(None,None,3),slice(None,None,3))
fig,ax=plt.subplots()
im=ax.imshow(phi,origin="lower",extent=[0,1,0,1],cmap='jet')
ax.quiver(X[skip],Y[skip],dx[skip],dy[skip])
fig.colorbar(im).set_label("$\phi/V$",rotation=270,labelpad=10)
plt.xlabel("$x/L$")
plt.ylabel("$y/L$")
ax.set(aspect=1)
plt.show()
#Gauss-seidel
def poisson_gauss(phi_c,phi_u,phi_d,phi_l,phi_r,roll,omega):
    return (1.+omega)*(phi_u+phi_d+phi_r+phi_l+roll)/4.-omega*phi_c

phi=np.zeros([101,101],dtype=np.float32)
phi[21:82,20]=1.
phi[21:82,80]=-1.
roll=np.zeros([101,101])
epi=1e-6
delta=epi+1.
omega=0.9
x=[i for i in range(100)]
x.append(-1)
x=np.array(x)
while delta>epi:
    i=j=1
    delta=phi[1,1]-poisson_gauss(phi[1,1],phi[1,2],phi[1,0],phi[0,1],phi[2,1],roll[1,1],omega)
    while i<100:
        while j<100:
            old=phi[i,j]
            if j==20 or j==80:
                if i>20 and i<80:
                    pass
                else:
                    phi[i,j]=poisson_gauss(phi[i,j],phi[i,j+1],phi[i,j-1],phi[i-1,j],phi[i+1,j],roll[i,j],omega)
            else:
                phi[i,j]=poisson_gauss(phi[i,j],phi[i,j+1],phi[i,j-1],phi[i-1,j],phi[i+1,j],roll[i,j],omega)
            if delta<np.abs(old-phi[i,j]):
                delta=np.abs(old-phi[i,j])
            j+=1
        j=1
        i+=1
    print(delta)
x=np.linspace(0,10,100)
y=np.linspace(0,10,100)
X,Y=np.meshgrid(x,y)
dy,dx=np.gradient(-phi)
skip=(slice(None,None,3),slice(None,None,3))
fig,ax=plt.subplots()
im=ax.imshow(phi,origin="lower",extent=[0,10,0,10],cmap='jet')
ax.quiver(X[skip],Y[skip],dx[skip],dy[skip])
fig.colorbar(im).set_label("$\phi/V$",rotation=270,labelpad=10)
plt.xlabel("$x/L$")
plt.ylabel("$y/L$")
ax.set(aspect=1)
plt.show()
