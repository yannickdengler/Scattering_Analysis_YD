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


if __name__ == "__main__":
    def generalized_mom(E_pipi, m_meson):
        if E_pipi < 2*m_meson:
            return 0
        return 2*np.arcsin(np.sqrt(0.5*(np.cosh(E_pipi/2.)-np.cosh(m_meson))))

    # # SASA PAPER ANALYSIS 2202.10110v4 (heavy ensemble)

    # m_D = 1927
    # m_D_star = 2049 
    # m_meson = (m_D+m_D_star)/2.

    # hbar_c = 197.32698

    # L_1 = 2.745/hbar_c
    # E_cm_1 = 1.02
    # P_1 = m_meson*generalized_mom(E_cm_1*2, 1)
    # q_1 = P_1 * L_1/(2*np.pi)
    # tan_PS_1 = np.pi**(3/2.)*q_1/Zeta(q_1**2)


    # L_2 = 2.048/hbar_c
    # E_cm_2 = 1.038
    # P_2 = m_meson*generalized_mom(E_cm_2*2, 1)
    # q_2 = P_2 * L_2/(2*np.pi)
    # tan_PS_2 = np.pi**(3/2.)*q_2/Zeta(q_2**2)




    # L_3 = 2.78/hbar_c
    # E_cm_3 = 1.022
    # P_3 = m_meson*generalized_mom(E_cm_3*2, 1)
    # q_3 = P_3 * L_3/(2*np.pi)
    # tan_PS_3 = np.pi**(3/2.)*q_3/Zeta(q_3**2)

    # print(P_1, P_2,P_3)
    # print(tan_PS_1, tan_PS_2,tan_PS_3)

    # plt.scatter((P_1/(2*m_meson))**2, P_1/(tan_PS_1*2*m_meson), label = "1")
    # plt.scatter((P_2/(2*m_meson))**2, P_2/(tan_PS_2*2*m_meson), label = "2")
    # plt.scatter((P_3/(2*m_meson))**2, P_3/(tan_PS_3*2*m_meson), label = "3")

    # plt.grid()
    # plt.legend()
    # plt.show()

    # exit()

    # SU(3) ANALYSIS

    # SET 1 b = 2.7984, k = 0.2984, g = 1.317, m_H = 148, a = 287
    # m_inf = 0.695162
    # E = [1.40474, 1.38779, 1.41455, 1.47866]
    # L = [12,14,10,8]
    # print(len(E), len(L), 2*m_inf)
    # P = []
    # P_2 = []
    # s = []
    # PS = []
    # q = []
    # Zeta_tmp = []
    # tan_PS = []
    # PS = []
    # E_plot = []

    # for i in range(len(E)):
    #     if E[i] < 2*m_inf:
    #         s.append(0)
    #         P.append(0)
    #         P_2.append(0)
    #         q.append(0)
    #         Zeta_tmp.append(0)
    #         tan_PS.append(0)
    #         PS.append(0)
    #     else:
    #         E_plot.append(E[i])
    #         s.append(E[i]*E[i])
    #         P.append(generalized_mom(E[i], m_inf))
    #         P_2.append(P[i]*P[i])
    #         q.append(P[i]*L[i]/(2*np.pi))
    #         Zeta_tmp.append(Zeta(q[i]**2))
    #         tan_PS.append(np.pi**(3/2.)*q[i]/Zeta_tmp[i])
    #         PS.append(np.arctan(tan_PS[i]))
    # print(E, tan_PS, P)

    # for i in range(len(E)):
    #     print(L[i], P[i], tan_PS[i])


    # plt.scatter(P, tan_PS, color = "red", label = "My data", s = 30)
    # plt.legend()
    # plt.show()
    # exit()
    # # BERND PAPER ANALYSIS

    # # SET 1 b = 2.7984, k = 0.2984, g = 1.317, m_H = 148, a = 287
    # lattice_constant = 287
    # m_inf = 0.280
    # E = [147,192, 236, 156, 183, 268, 150, 177, 149, 173, 148, 163, 183, 193, 257, 147, 171, 295, 149, 164, 219]
    # for i in range(len(E)):
    #     E[i] = E[i]/lattice_constant
    # L = [8,8,8,12,12,12,16,16,20,20,24,24,24,24,24,28,28,28,32,32,32]

    # # with open("/home/dengler_yannick/Documents/Wolfram Mathematica/Axel_Scattering/energy_level_data_Bernd_short.dat", "w") as file:
    # #     for i in range(len(E)):
    # #         file.write("IDC\t\t")
    # #         file.write("%1.3f\t"%1)
    # #         file.write("%1.3f\t"%1)
    # #         file.write("%1.3f\t"%1)
    # #         file.write("%i\t"%1)
    # #         file.write("%i\t"%L[i])
    # #         file.write("%i\t"%0)
    # #         file.write("%e\t"%0)
    # #         file.write("%e\t"%m_inf)
    # #         file.write("%e\t"%0)
    # #         file.write("%e\t"%E[i])
    # #         file.write("%e\t"%0)
    # #         file.write("\n")
    # # exit()
    # print(len(E), len(L), 2*m_inf)
    # P = []
    # P_2 = []
    # s = []
    # PS = []
    # q = []
    # Zeta_tmp = []
    # tan_PS = []
    # PS = []
    # E_plot = []

    # for i in range(len(E)):
    #     # if E[i] < 2*m_inf:
    #     #     s.append(0)
    #     #     P.append(0)
    #     #     P_2.append(0)
    #     #     q.append(0)
    #     #     Zeta_tmp.append(0)
    #     #     tan_PS.append(0)
    #     #     PS.append(0)
    #     # else:
    #     E_plot.append(E[i])
    #     s.append(E[i]*E[i])
    #     P.append(generalized_mom(E[i], m_inf))
    #     P_2.append(P[i]*P[i])
    #     q.append(P[i]*L[i]/(2*np.pi))
    #     Zeta_tmp.append(Zeta(q[i]**2))
    #     tan_PS.append(np.pi**(3/2.)*q[i]/Zeta_tmp[i])
    #     PS.append(np.arctan(tan_PS[i]))
    # print(E, tan_PS, P)

    # for i in range(len(E)):
    #     print(P[i], PS[i])


    # E_bernd = [172.711,177.113,183.451,192.254,218.662,235.915,256.69,267.606,294.894]
    # tan_PS_bernd = [-0.457,-0.343,-0.257,-0.114,0.629,-0.543,0.457,6.229,0.943]

    # for i in range(len(E)):
    #     E[i] = E[i]*lattice_constant
    # plt.scatter(E, tan_PS, color = "red", label = "My data", s = 30)
    # plt.scatter(E_bernd, tan_PS_bernd, color = "blue", label = "Bernd", s = 20)
    # plt.legend()
    # plt.show()





    # exit()




















    # m_meson = 0.1673
    # E = 0.5107
    # P = generalized_mom(E, m_meson)
    # print(P)
    # s = E*E
    # print(s)
    # # P = 0.193964
    # q = P*16/(2*np.pi)
    # Zeta_tmp = Zeta(q**2)
    # tan_PS = np.pi**(3/2.)*q/Zeta_tmp
    # print(tan_PS)
    # print(np.arctan(tan_PS))
    # print(np.arctan(tan_PS)*360/(2*np.pi))
    # print(P**3/(np.sqrt(s)*tan_PS))




    # print(Zeta(q=-0.095900719))
    # print(Zeta(q=0.472894247))
    # print(Zeta(q=1.441591313))
    # print(Zeta(q=2.627007612))
    # # for q in range(10):
    # #     print(q, Zeta(q))
    # # exit()

    xarr1 = []
    yarr1 = []
    xarr2 = []
    yarr2 = []
    xarr3 = []
    yarr3 = []
    xarr4 = []
    yarr4 = []
    yarr41 = []

    orders = ((7,20,4), (8,21,5), (9,22,6), (10,23,7))

    for q in np.linspace(-1, 1.5, 1000):
    # for q in np.linspace(-1, 10.2, 400):
    # for q in np.linspace(-0.1, 1.1, 400):
        Zeta_val = Zeta(q)
        print(q, Zeta_val)
        # xarr1.append(q)
        # yarr1.append(func(q, orders[0]))
        # xarr2.append(q)
        # yarr2.append(func(q, orders[1]))
        # xarr3.append(q)
        # yarr3.append(func(q, orders[2]))
        xarr4.append(q)
        yarr4.append(Zeta_val)
        yarr41.append(q/Zeta_val)

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


    plt.plot(xarr4, yarr4, label = "Zeta")
    plt.plot(xarr4, yarr41, label = "q/Zeta")

    plt.legend()

    plt.savefig("Zeta.pdf")

    plt.show()                                                  # Function is not independant on order. Check for the error somewhere!!!

