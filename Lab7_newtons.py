import numpy as np
import matplotlib.pylab as plt
import Lab6_nonlinear_relaxation as relax

def binary_search(f,interval=[0,10],i=0,epi=10e-7):
    x=float(interval[0])+interval[1]
    x/=2
    df=f(x+epi)-f(x-epi)
    df/=epi
    if f(x)+epi<0:
        if df>0:
            return binary_search(f,[x,interval[1]],i+1,epi)
        elif df<0:
            return binary_search(f,[interval[0],x],i+1,epi)
        else:
            df=f(x+2*epi)-f(x)
            df/=epi
            if df>0:
                return binary_search(f,[x,interval[1]],i+1,epi)
            elif df<0:
                return binary_search(f,[interval[0],x],i+1,epi)
            else:
                return False
    elif f(x)-epi>0:
            if df<0:
                return binary_search(f,[x,interval[1]],i+1,epi)
            elif df>0:
                return binary_search(f,[interval[0],x],i+1,epi)
            else:
                df=f(x+2*epi)-f(x)
                df/=epi
                if df<0:
                    return binary_search(f,[x,interval[1]],i+1,epi)
                elif df>0:
                    return binary_search(f,[interval[0],x],i+1,epi)
                else:
                    return False
    elif np.abs(f(x))<epi:
        return x,i
def newton(f,epi,degree,interval):
    x=interval[0]
    df=f(x+1e-5)-f(x-1e-5)
    df/=(2e-5)
    d2f=f(x+1e-5)-2*f(x)+f(x-1e-5)
    d2f/=1e-10
    tmp=[x,x-f(x)/df]
    err=np.abs(tmp[1]-tmp[0])
    step=[]
    ans=[]
    no_of_root=0
    while no_of_root<degree :
        step.append(0)
        step[-1]+=1
        err=np.abs(tmp[1]-tmp[0])
        print(err)
        while err>epi:
            tmp[0]=tmp[1]
            df=f(tmp[0]+1e-5)-f(tmp[0]-1e-5)
            df/=(2e-5)
            d2f=f(tmp[0]+epi)-2*f(tmp[0])+f(tmp[0]-epi)
            d2f/=1e-10
            try:
                print(df,d2f)
            except ValueError:
                return False
            tmp[1]=tmp[0]-f(tmp[0])/df
            step[-1]+=1
            err=np.abs(((-1.*d2f)/(2.*df)))*err**2
            if np.abs(tmp[1])>1e200:
                ans.append("NA")
        if len(ans)>no_of_root:
            pass
        else:
            ans.append(tmp[1])
        no_of_root+=1
        x+=(interval[0]+interval[1])/float(degree)
        tmp[0]=x
        df=f(tmp[0]+epi)-f(tmp[0]-epi)
        df/=(2*epi)
        if df!=0:
            tmp[1]=tmp[0]-f(tmp[0])/df
        else:
            x+=epi
            tmp[0]=x
            df=f(tmp[0]+epi)-f(tmp[0]-epi)
            df/=epi
            tmp[1]=tmp[0]-f(tmp[0])/df
    return ans,step
def f_1(x):
    return 5*np.math.exp(-1*x)+x-5
def f_1_relax(x):
    return 5-5*np.exp(-1*x)
def polynomial(x):
    return 924*x**6-2772*x**5+3150*x**4-1680*x**3+420*x**2-42*x+1
'''
root,i=binary_search(f_1,[0,10])
print("Binary search method: ",root,i)
root,i=relax.relaxation(f_1_relax,10,10e-7)
print("Relaxation",root,i)
ans,step=newton(f_1,10e-7,2,[0,10])
print("newton's method:",ans,step)
x_axis=np.linspace(0,1)
y=[polynomial(x) for x in x_axis]
plt.plot(x_axis,y)
plt.show()
ans,step=newton(polynomial,10e-12,6,[0,1])
print("newton's method:",ans,step)
ans,step=binary_search(polynomial,[0,1],0,10e-12)
print("Binary search:",ans,step)
'''
