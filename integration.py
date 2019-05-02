import numpy as np
import gaussxw
def trapezoidal(f,a,b,no_of_step):
    a=float(a)
    b=float(b)
    step=(b-a)/no_of_step
    summation=0
    summation+=(f(a)+f(b))*0.5
    for i in range(1,no_of_step):
            summation+=f(a+i*step)
    return summation*step

def simpson(f,a,b,no_of_step):
    a=float(a)
    b=float(b)
    step=(b-a)/no_of_step
    summation=0
    summation+=(f(a)+f(b))
    for i in range(0,no_of_step):
        summation+=4*(f((2*a+(2*i+1)*step)/2))
    for i in range(1,no_of_step):
        summation+=2*(f(a+i*step))
    return summation*(step/6)

def gauss_quad(gauss_x,gauss_w,lower_bound,upper_bound):
    return 0.5*(upper_bound-lower_bound)*gauss_x+0.5*(lower_bound+upper_bound),0.5*(upper_bound-lower_bound)*gauss_w

class numerical_method():
    def __init__(self,method,adaptive_method,f,error_evaluate,a,b,step,error=10e-7):
        self.res=[]
        self.a=float(a)
        self.b=float(b)
        self.step=step
        self.error=error
        self.f=f
        self.h=(self.b-self.a)/self.step
        self.method=method
        self.adaptive_method=adaptive_method
        self.adaptive_res=[]
        self.adaptive_step=0
        self.no_of_sample=1
        self.res.append(self.method(self.f,self.a,self.b,self.step))
        self.step*=2
        self.res.append(self.method(self.f,self.a,self.b,self.step))
        self.error_evaluate=error_evaluate
        self.accuracy=np.abs(self.error_evaluate(self.res[-1],self.res[-2]))
    def increase_accuracy(self):
        if len(self.res)>50:
            self.step*=2
        self.res.append(self.method(self.f,self.a,self.b,self.step))
        self.accuracy=self.error_evaluate(self.res[-1],self.res[-2])
    def get_accuracy(self):
        self.accuracy=np.abs(self.error_evaluate(self.res[-1],self.res[-2]))
        return self.accuracy
    def get_adaptive_accuracy(self):
        return np.abs(self.error_evaluate(self.adaptive_res[-1],self.adaptive_res[-2]))
    def change_interval(self,coord,step=10):
        self.res=[]
        self.step=10
        self.x,self.y=coord[0],coord[1]
        self.res.append(self.method(self.f,self.a,self.b,self.step))
        self.increase_accuracy()
        self.accuracy=self.get_accuracy()
    def latest_res(self):
        return self.res[-1]
    def adaptive_cal(self,a,b,error=10e-7,mode=0):
        if(self.adaptive_method==None):
            self.adaptive_res=[]
            self.adaptive_res.append(self.method(self.f,a,b,int(self.step/2)))
            self.adaptive_res.append(self.method(self.f,a,b,self.step))
            self.accuracy=self.get_adaptive_accuracy()
            while(self.accuracy)>=error:
                self.step*=2
                self.h=(b-a)/self.step
                extra_term=self.h*np.sum([self.f(a+x*self.h) for x in range(1,self.step,2)])
                self.adaptive_res.append(0.5*self.adaptive_res[-1]+extra_term)
                self.accuracy=self.get_adaptive_accuracy()
            self.res.append(self.adaptive_res[-1])
            if(mode==1):
                sum=0
                graph=[]
                x=[]
                h=(b-a)/self.step
                for item in range(0,self.step):
                    x.append(a+item*h)
                    sum+=(self.f(a+item*h)+self.f(a+(item+1)*h))*(h/2)
                    graph.append(sum)
                return (x,graph)
        if(self.adaptive_method!=None):
                c=(a+b)/2.0
                i=0
                sample={0:[a,c,b]}
                self.no_of_sample=3
                self.adaptive_res=[0]
                tmp_sum=0
                tmp_last=0
                for coord in range(1,len(sample[i])-1):
                    left=self.method(self.f,sample[i][coord-1],sample[i][coord],1)
                    right=self.method(self.f,sample[i][coord],sample[i][coord+1],1)
                    whole=self.method(self.f,sample[i][coord-1],sample[i][coord+1],1)
                    tmp_sum+=(left+right)
                    tmp_last+=whole
                self.adaptive_res.append(tmp_sum)
                self.accuracy=(tmp_sum-tmp_last)/15.0
                while (self.accuracy)>=error:
                    i+=1
                    self.no_of_sample+=2**i
                    new_sample=[]
                    sample[i]=[]
                    tmp_sum=0
                    tmp_last=0
                    for item in range(1,len(sample[i-1])-1,2):
                        new_sample.append((sample[i-1][item-1]+sample[i-1][item])/2.0)
                        new_sample.append((sample[i-1][item]+sample[i-1][item+1])/2.0)
                    sample[i].append(sample[i-1][0])
                    for coord in range(0,len(new_sample)):
                        sample[i].append(new_sample[coord])
                        sample[i].append(sample[i-1][coord+1])
                    for coord in range(0,len(sample[i-1])-1):
                        left=self.method(self.f,sample[i-1][coord],new_sample[coord],1)
                        right=self.method(self.f,new_sample[coord],sample[i-1][coord+1],1)
                        tmp_sum+=(left+right)

                    self.adaptive_res.append(tmp_sum)
                    self.accuracy=np.abs(tmp_sum-self.adaptive_res[-2])/15.0
                self.res.append(self.adaptive_res[-1])
