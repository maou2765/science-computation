def simpsons_rule(f,a,b):
    c = (a+b) / 2.000
    h3 = abs(b-a) / 6.0
    return h3*(f(a) + 4.0*f(c) + f(b))

def recursive_asr(f,a,b,eps,whole):
    "Recursive implementation of adaptive Simpson's rule."
    c = (a+b) / 2.0
    left = simpsons_rule(f,a,c)
    right = simpsons_rule(f,c,b)
    print("c=",c," a=",a," b=",b," left=",left," right=",right," whole=",whole)
    if abs(left + right - whole) <= 15*eps:
        print(left + right + (left + right - whole)/15.0)
        return left + right + (left + right - whole)/15.0
    return recursive_asr(f,a,c,eps/2.0,left) + recursive_asr(f,c,b,eps/2.0,right)

def adaptive_simpsons_rule(f,a,b,eps):
    "Calculate integral of f from a to b with max error of eps."
    return recursive_asr(f,a,b,eps,simpsons_rule(f,a,b))

from math import sin
print (adaptive_simpsons_rule(sin,0,2,.000000001))
