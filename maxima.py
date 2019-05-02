import numpy as np
def f(x):
    return (np.math.sin(x)/x)**2
err=1.
epi=0.01
upd=5.
lob=4.
i=0
while err>epi:
    i+=1
    mid=(upd+lob)/2
    fup=f((upd+mid)/2.)
    flo=f((lob+mid)/2.)
    if fup>flo:
        lob=(lob+mid)/2.
    else:
        upd=(upd+mid)/2.
    err=upd-lob
    print(fup,flo,upd,lob)
print(fup,mid)
print(i)
i=0
x=[4.,4.25,4.5,5.]
x[2]=(x[3]-x[0])/1.618+x[0]
x[1]=x[3]-x[2]+x[0]
err=x[2]-x[1]
while err>epi:
    i+=1
    if f(x[1])>f(x[2]):
        x[3]=x[2]
        x[2]=(x[3]-x[0])/1.618+x[0]
        x[1]=x[3]-x[2]+x[0]
    elif f(x[2])>f(x[1]):
        x[0]=x[1]
        x[2]=(x[3]-x[0])/1.618+x[0]
        x[1]=x[3]-x[2]+x[0]
    err=x[2]-x[1]
    print(x)
print(f(x[2]),x[2])
print(i)
