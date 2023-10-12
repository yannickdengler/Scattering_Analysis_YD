import numpy as np
import scipy.integrate as scpint
import matplotlib.pyplot as plt
from decimal import Decimal
from decimal import getcontext

# getcontext().prec = 40

def y_elem_Z3(order):
    y = []
    for i in [x for x in range(-order, order+1)]:
        for j in [x for x in range(-order, order+1) if abs(x)+abs(i) < order]:
            for k in [x for x in range(-order, order+1) if abs(x)+abs(i)+abs(j) < order]:
                y.append((i,j,k))
    return y

def y_squared(order):
    res = []
    for y in y_elem_Z3(order):
        res.append(y[0]*y[0]+y[1]*y[1]+y[2]*y[2])
    return res
    
def term1(q, order=10):
    res = 0
    for y_2 in y_squared(order):
        res += np.exp(-(y_2-q))/(y_2-q)
    return res

def factorial(x):
    res = 1
    for i in range(1, x+1):
        res = res*i
    return res

def term2(q, order=23):
    res = 0
    for k in range(order):
        res += (q**k)/(factorial(k)*(k-0.5))
    return np.pi**(3/2.)*res


def exp_3(u, order):
    res = 0
    for y_2 in y_squared(order):
        res += np.exp(-np.pi**2*y_2/u)
    return res-1

def e_u_q2(u, q, order):
    return np.exp(u*q)*(np.pi/u)**(3/2.)*exp_3(u, order)

def term3(q, order=7):
    return scpint.quad(func=e_u_q2, a=0, b=1, args=(q, order,))[0]


def Zeta(q, orders=(10,23,7)):
    return (term1(q, orders[0]) + term2(q, orders[1]) + term3(q, orders[2]))/(np.sqrt(4*np.pi))

# print(Zeta(q=-0.095900719))
# print(Zeta(q=0.472894247))
# print(Zeta(q=1.441591313))
# print(Zeta(q=2.627007612))
# # for q in range(10):
# #     print(q, Zeta(q))
# # exit()

# xarr1 = []
# yarr1 = []
# xarr2 = []
# yarr2 = []
# xarr3 = []
# yarr3 = []
# xarr4 = []
# yarr4 = []

# # orders = ((7,20,4), (8,21,5), (9,22,6), (10,23,7))

# for q in np.linspace(-1, 10.2, 400):
# # for q in np.linspace(-0.1, 1.1, 400):
#     print(q)
#     func = Zeta
# #     # xarr1.append(q)
# #     # yarr1.append(func(q, orders[0]))
# #     # xarr2.append(q)
# #     # yarr2.append(func(q, orders[1]))
# #     # xarr3.append(q)
# #     # yarr3.append(func(q, orders[2]))
#     xarr4.append(q)
#     yarr4.append(func(q))

# # # # orders = (20,21,22,23)

# # # # for q in np.linspace(-1, 10.2, 10000):
# # # #     print(q)
# # # #     func = term2
# # # #     xarr1.append(q)
# # # #     yarr1.append(func(q, orders[0]))
# # # #     xarr2.append(q)
# # # #     yarr2.append(func(q, orders[1]))
# # # #     xarr3.append(q)
# # # #     yarr3.append(func(q, orders[2]))
# # # #     xarr4.append(q)
# # # #     yarr4.append(func(q, orders[3]))

# # # # plt.yscale("log")

# plt.ylim((-10,10))

# # plt.plot(xarr1, yarr1, label = str(orders[0]))
# # plt.plot(xarr2, yarr2, label = str(orders[1]))
# # plt.plot(xarr3, yarr3, label = str(orders[2]))
# # plt.plot(xarr4, yarr4, label = str(orders[3]))


# plt.plot(xarr4, yarr4, label = "Zeta")

# plt.legend()

# plt.savefig("Zeta.pdf")

# plt.show()                                                  # Function is not independant on order. Check for the error somewhere!!!

