# -*- coding: utf-8 -*-
from sympy import *
init_printing(use_unicode=False, wrap_line=False, no_global=True)
x = Symbol('x')
sigma=1
mu=1
f=1/sigma/sqrt(2*pi)*exp(-1/2*(x-mu)**2/sigma**2)
integ=integrate(f, (x,-1,2))
print(float(integ))

