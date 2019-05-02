import numpy as np
import matplotlib.pyplot as plt
import datetime
My_time_start=datetime.datetime.now()
def voltage(coord):#calculate the Madelung constant
    distance=0
    oddoreven=0
    distance=((coord[0])**2.0+(coord[1])**2.0+(coord[2])**2.0)**(-0.5)
    oddoreven=(coord[0]+coord[1]+coord[2])
    oddoreven=oddoreven%2
    if distance==0:
        return 0
    elif oddoreven:
        return (distance*(-1))
    else:
        return distance

def degeneracy(Elevel):#this function is used to calculate degenercy of each voltage level
    if Elevel[0]==Elevel[1]==Elevel[2]:
        return 8 #2^3
    elif Elevel[1]==Elevel[2]==0:
        return 6#3*2^1
    elif Elevel[0]==Elevel[1] and Elevel[2]==0:
        return 12#3C2*2^2
    elif Elevel[0]!=Elevel[1] and Elevel[2]==0:
        return 24#3*2!*2^2
    elif Elevel[2]!=0 and Elevel[0]==Elevel[1]:
        return 24#3*2^3
    elif Elevel[1]==Elevel[2] and Elevel[2]!=0:
        return 24#3*2^3
    elif Elevel[0]==Elevel[2] and Elevel[1]!=0:
        return 24#3*2^3
    elif Elevel[1]!=Elevel[2] and Elevel[0] !=Elevel[1] and Elevel[0] != Elevel[2] and Elevel[2]!=0 and Elevel[1]!=0:
        return 48#3!*2^3

#coord_that_need2cal=[] #create a list consist of coordinate corresipond to each L
#Lcoord_that_need2cal={}#a dict that contain the above list
element_in_the_list=[0,0,0]#contain coordinate
L=range(1,251)
Madelung=[]
My_time_start=datetime.datetime.now()
voltage_total=0
for x in range(L[-1],0,-1):
    element_in_the_list[0]=x
    if element_in_the_list[1]==0:
        for y in range(x,-1,-1):
            element_in_the_list[1]=y
            if element_in_the_list[2]==0:
                for z in range(y,-1,-1):
                    element_in_the_list[2]=z
                    mul=degeneracy(element_in_the_list)
                    voltage_total+=(mul*voltage(element_in_the_list))#Add all voltage and corr. degenercy to get correct Madelung of given L
    Madelung.append(voltage_total)
                    #coord_that_need2cal.append(element_in_the_list.copy())
    #Lcoord_that_need2cal[x]=coord_that_need2cal.copy()
    #coord_that_need2cal=[]
Your_time_start=datetime.datetime.now()
My_time_need_list=Your_time_start-My_time_start
print("time elapsed for generation of the list: ", My_time_need_list)
print('L=',499,'M=',Madelung[-2],' L=',500,'M=',Madelung[-1])
'''
largest_L=L[0]
for term in L:
    if term >largest_L:
        largest_L=term
voltage_total=0
My_time_start=datetime.datetime.now()
for j in range(1,largest_L+1):
        for coord in Lcoord_that_need2cal[j]:
            mul=degeneracy(coord)
            voltage_total+=(mul*voltage(coord))#Add all voltage and corr. degenercy to get correct Madelung of given L
        Madelung.append(voltage_total)
print('L= ',largest_L,'M',Madelung[-1])
Your_time_start=datetime.datetime.now()
My_time_need=Your_time_start-My_time_start
print("time elapsed: ", My_time_need+My_time_need_list)
log_M=[np.log(np.negative(M)) for M in Madelung]
log_L=[np.log(L) for L in L]

plt.plot(log_L, log_M,'bo-')
plt.xlabel('log(L)')
plt.ylabel('log(M)')
plt.show()
'''
