import numpy as np
import matplotlib.pylab as plt
def my_cos_1(x):
    return np.math.cos(x-20*(np.math.sin(x)))
def my_linear(a,b,f):
    slope=(f(b)-f(a))/float(b-a)
    def line(x):
        return (slope*x+(f(a)-slope*a))
    return line
x=[(np.math.pi/1000)*x for x in range(0,1000)]
my_cos=[my_cos_1(x_coord) for x_coord in x]
plt.plot(x,my_cos)
start=end=0
for n in range(1,21):
    end=start+20*np.math.pi/400
    line_1=my_linear(start,end,my_cos_1)
    x=[(np.math.pi/400)*x+start for x in range(0,21)]
    my_line=[line_1(x_coord) for x_coord in x]
    print(end)
    start=end
    plt.plot(x,my_line)
plt.show()
