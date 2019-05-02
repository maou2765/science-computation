import numpy as np
import integration as ig
import gaussxw as gauss
import matplotlib.pylab as plt
def F_1(x):
    return np.reciprocal(np.sqrt(np.math.sin(x)))
No_Sample_list=[10,50,100,1000,5000]
res_F_1=[]
for N in No_Sample_list:
    gauss_x,gauss_w=gauss.gaussxw(N)
    gauss_x_N,gauss_w_N=ig.gauss_quad(gauss_x,gauss_w,0,1)
    sum=0
    for i in range(N):
        sum+=F_1(gauss_x_N[i])*gauss_w_N[i]
    res_F_1.append(sum)
for i in range(len(res_F_1)):
    print("When the number of sample is ",No_Sample_list[i],', the result is',res_F_1[i])
F_1_trap=ig.numerical_method(ig.trapezoidal,None,F_1,lambda x,y:(x-y)/3.0,0.0001,1,10,10e-7)

for x in [0.2,0.1,0.01,0.001,0.00000001]:
    F_1_trap.adaptive_cal(x,1,10e-7,0)
    print("I(",x,")=",F_1_trap.res[-1])
    
'''
x_coord,y_coord=F_1_trap.adaptive_cal(0.001,1,10e-7,1)
y_coord=[y_coord[-1]-y for y in y_coord]
x_coord=[np.math.log(x) for x in x_coord]
plt.plot(x_coord,y_coord)
plt.show()
'''
def F_2(x):
    if x==0:
        return 0
    return 1/(np.sqrt(np.math.sin(x)))-1/(np.sqrt(x))

def F_3(x):
    return np.reciprocal(np.sqrt(x))
F_2_trap=ig.numerical_method(ig.trapezoidal,None,F_2,lambda x,y:(x-y)/3.0,0,1,10,10e-7)
res_F_2_trap=[]
F_2_trap.adaptive_cal(0,1,10e-8,0)
print("I0=",F_2_trap.res[-1],"with error: ",F_2_trap.accuracy)
print("I=",F_2_trap.res[-1]+2)
def F(x,y,kbT):
    return np.negative(kbT*np.math.log(1+np.math.exp(-2*(np.math.cos(x)+np.math.cos(y))/kbT)))
gauss_x,gauss_w=gauss.gaussxw(100)
guass_x,gauss_w=ig.gauss_quad(gauss_x,gauss_w,np.negative(np.math.pi),np.math.pi)
sample_pt=[]
weight=[]
for x in gauss_x:
    for y in gauss_x:
        sample_pt.append([x,y])
for w1 in gauss_w:
    for w2 in gauss_w:
        weight.append([w1,w2])
def Free_energy(sample_pt,weight,kbT):
    sum=0
    for i in range(len(sample_pt)):
        sum+=(F(sample_pt[i][0],sample_pt[i][1],kbT)*weight[i][0]*weight[i][1])
    return sum
kbT_list=[x*0.1 for x in range(1,101)]
kbt_1=Free_energy(sample_pt,weight,1)
print("Free energy when KbT=1: ",kbt_1)
free_energy=[Free_energy(sample_pt,weight,x) for x in kbT_list]
plt.plot(kbT_list,free_energy)
plt.show()
