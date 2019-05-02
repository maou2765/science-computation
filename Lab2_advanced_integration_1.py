import integration as ig
import numpy as np
import matplotlib.pylab as plt
import gaussxw
def f_1(t):
    if t==0:
        return 0
    else:
        return ((np.math.sin(10*t))**2.0)/t

error=10e-7
no_of_step=10
x=2
E_trap=ig.numerical_method(ig.trapezoidal,None,f_1,lambda x,y: (x-y)/3.0,0,x,no_of_step,error=error)
E_simpson=ig.numerical_method(ig.simpson,1,f_1,lambda x,y: (x-y)/15.0,0,x,no_of_step,error)
E_trap.adaptive_cal(0,E_trap.b)
print("Adaptive trapezoidal rule:",E_trap.latest_res(),"now accuracy",E_trap.get_adaptive_accuracy(),"with step: ",E_trap.step)
E_simpson.adaptive_cal(0,E_simpson.b,10e-7)
print("Adaptive Simpson's rule:",E_simpson.latest_res(),"now accuracy",E_simpson.get_adaptive_accuracy(),"with step: ",E_simpson.no_of_sample)

x_coord,y_coord=E_trap.adaptive_cal(0,3,10e-7,1)
#x_coord=[x*E_trap.h for x in range(0,E_trap.step+1)]
plt.plot(x_coord,y_coord)
plt.show()
def C_V(T):
    def my_CV(x):
        if x==0:
            return 0
        else:
            answer=1
            answer*=9*(10**(-3.0))*6.022*10**28*1.38*10**(-23)
            answer*=np.power(T/428.0,3)
            answer*=np.product([x**4,np.exp(x)])
            answer*=np.reciprocal(np.square(np.exp(x)-1))
            return answer
    return my_CV
C_V_500=C_V(500)
C_V_500_trap=ig.numerical_method(ig.trapezoidal,None,C_V_500,lambda x,y:(x-y)/3.0,0,428.0/500,no_of_step,error=error)
C_V_500_trap.adaptive_cal(0,C_V_500_trap.b)
print("CV at 500k","Adaptive trapezoidal rule:",C_V_500_trap.latest_res(),"now accuracy",C_V_500_trap.get_adaptive_accuracy(),"with step: ",C_V_500_trap.step)
C_V_500_simpson=ig.numerical_method(ig.simpson,1,C_V_500,lambda x,y:(x-y)/3.0,0,428.0/500,no_of_step,error=error)
C_V_500_simpson.adaptive_cal(0,C_V_500_simpson.b)
print("CV at 500k","Adaptive Simpson's rule:",C_V_500_simpson.latest_res(),"now accuracy",C_V_500_simpson.get_adaptive_accuracy(),"with step: ",C_V_500_simpson.no_of_sample)

C_V_5=C_V(5)
C_V_5_trap=ig.numerical_method(ig.trapezoidal,None,C_V_5,lambda x,y:(x-y)/3.0,0,428.0/5,no_of_step,error=error)
C_V_5_trap.adaptive_cal(0,C_V_5_trap.b)
print("CV at 5k","Adaptive trapezoidal rule:",C_V_5_trap.latest_res(),"now accuracy",C_V_5_trap.get_adaptive_accuracy(),"with step: ",C_V_5_trap.step)
C_V_5_simpson=ig.numerical_method(ig.simpson,None,C_V_5,lambda x,y:(x-y)/15.0,0,428.0/5,no_of_step,error=error)
C_V_5_simpson.adaptive_cal(0,C_V_5_simpson.b)
gauss_x,gauss_w=gaussxw.gaussxw(50)
def CV_gauss(gauss_x,gauss_w,T,f):
    x,w=ig.gauss_quad(gauss_x,gauss_w,0,428.0/T)
    sum=0
    for i in range(50):
        sum+=f(x[i])*w[i]
    return sum
CV_500_gauss_value=CV_gauss(gauss_x,gauss_w,500,C_V_500)
CV_5_gauss_value=CV_gauss(gauss_x,gauss_w,5,C_V_5)
print("CV at 500k, Gaussian quadrature for 50 samples: ",CV_500_gauss_value,"which has ",(CV_500_gauss_value-C_V_500_trap.latest_res())/C_V_500_trap.latest_res()," fractional difference with the the result calculated by trapezoidal rule")
print("CV at 5k, Gaussian quadrature for 50 samples: ",CV_5_gauss_value,"which has ",(CV_5_gauss_value-C_V_5_trap.latest_res())/C_V_5_trap.latest_res()," fractional difference with the the result calculated by trapezoidal rule")


x_coord=[5+x*0.1 for x in range(4951)]
y_coord=[CV_gauss(gauss_x,gauss_w,x,C_V(x)) for x in x_coord]
plt.plot(x_coord,y_coord)
plt.show()
