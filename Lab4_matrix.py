import numpy.linalg as npla
import numpy as np

def SL_Vi(i):
    if i>=4:
        A=' 3. -1 -1'+(i-3)*' 0'+';'
        A+=' -1. 4 -1 -1'+(i-4)*' 0'+';'
        repeat_part=' -1. -1 4 -1 -1'
        j=0
        while True:
            if j<=(i-5):
                A+=j*' 0'+repeat_part+' 0'*(i-j-5)+';'
                j+=1
            else:
                break
        A+=' 0'*(i-4)+' -1 -1 4 -1.;'
        A+=' 0'*(i-3)+' -1 -1 3.'
    elif i==3:
        A=' 3. -1 -1;'
        A+=' -1. 4 -1;'
        A+=' -1. -1 3'
    elif i<3:
        return False
    return np.array(np.mat(A),subok=True)
def b_generation(i):
    b='5. 5'
    b+=(i-2)*' 0'
    return np.array(np.mat(b),subok=True).T
def show_ans(i):
    A=SL_Vi(i)
    b=b_generation(i)
    x=npla.solve(A,b)
    j=0
    print('N=',i)
    for V in x:
        j+=1
        print('V'+str(j),'=',V)
show_ans(3)
show_ans(6)
def band_solver(A,b,no_of_diagonal=0):
    col=0
    row=0
    if A.shape[0]!=A.shape[1] or A.shape[1]!=b.shape[0]:
        return False
    if no_of_diagonal==0:
        res=npla.solve(A,b)
        return res
    while col<A.shape[1]:
            div=1/A[row,col]
            A[row]*=div
            b[row]*=div
            if (row+no_of_diagonal)>=A.shape[0]:
                for j in range(row+1,A.shape[0]):
                    tmp_mul=A[j,col]/A[row,col]
                    for k in range((row+no_of_diagonal)-A.shape[0]):
                        A[j,col+k]-=(A[row,col+k]*tmp_mul)
                        b[j]-=(b[row]*tmp_mul)
            else:
                for j in range(row+1,row+no_of_diagonal):
                    tmp_mul=A[j,col]/A[row,col]
                    for k in range(no_of_diagonal):
                        A[j,col+k]-=(A[row,col+k]*tmp_mul)
                        b[j]-=(b[row]*tmp_mul)
            row+=1
            col+=1
    res=[b[-1]]
    for n in range(b.shape[0]-2,-1,-1):
        tmp=b[n]
        if no_of_diagonal==(A.shape[0]+5):
            for j in range(A.shape[1]-1,col-2,-1):
                tmp-=A[n,j]*res[A.shape[1]-1-j]
                tmp/=A[n,col-2]
        else:
            if (A.shape[1]-1-(col-2))<=no_of_diagonal:
                for j in range(A.shape[1]-1,col-2,-1):
                    tmp-=A[n,j]*res[A.shape[1]-1-j]
                    tmp/=A[n,col-2]
            else:
                for j in range(col-2+no_of_diagonal,col-2,-1):
                    tmp-=A[n,j]*res[A.shape[1]-1-j]
                    tmp/=A[n,col-2]
        col-=1
        res.append(tmp)
    return res
V6=band_solver(SL_Vi(6),b_generation(6),3)
V3=band_solver(SL_Vi(3),b_generation(3),3)
print('Below is solve by band_solver')
j=6
for item in V6:
    print('V'+str(j)+' =',item)
    j-=1
j=3
for item in V3:
        print('V'+str(j)+' =',item)
        j-=1
A=np.array(np.matrix('1. 2 3; 4 10 6; 7 8 9'),subok=True)
b=np.array(np.matrix('10. 11 12'),subok=True).T
test_band_solver=band_solver(A.copy(),b.copy())
j=3
for item in test_band_solver:
        print('V'+str(j)+' =',item)
        j-=1
print("solved by numpy")
test_standard=npla.solve(A.copy(),b.copy())
j=3
for item in test_standard:
        print('V'+str(j)+' =',item)
        j-=1


V10000=band_solver(SL_Vi(10000),b_generation(10000),3)
j=10000
for item in V10000:
    print('V'+str(j)+' =',item)
    j-=1

#14 using standard method
#6mins needed for N=10000
