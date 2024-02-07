import matplotlib.pyplot as plt
import error_classes as errcl 
import numpy as np
import scipy.integrate as scpint
from scipy.interpolate import interp1d as inter
import math
import h5py
import os

lightspeed = 299792.458
color_arr = ["blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick","blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick","blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick",]

def create_plot_file(name, data):
    with h5py.File("/home/dengler_yannick/Documents/Scattering_Analysis_YD/output/plot_files/"+name+".hdf5","w") as f:
        for key, val in data.items():
            f.create_dataset(key, data = val)

def read_plot_file(name):
    data = {}
    filename = "/home/dengler_yannick/Documents/Scattering_Analysis_YD/output/plot_files/"+name+".hdf5"
    if os.path.exists(filename):
        with h5py.File(filename,"r") as f:
            for key in f.keys():
                # print(type(f[key]))
                data[key] = f[key][()]
        return data



def P_cot_delta(P, coeffs):                                 # takes and return real unitless P = P/mass_Goldstone
    res = 0
    for i in range(len(coeffs)):
        if i == 0:
            res += coeffs[i]
        else:
            res += coeffs[i]*P**(2*i)
    return res

def get_data_from_file(file, order=2):
    coeff_strs = ["a", "b", "c", "d", "e"]
    fit_phase_shift = errcl.measurement(file)
    fit_phase_shift.read_from_HDF()
    P_cot_data = []
    coeffs = []

    # print(fit_phase_shift.infos.keys())
    beta = fit_phase_shift.infos["beta"]
    N_Ls = fit_phase_shift.infos["N_Ls"]
    m = fit_phase_shift.infos["m_1"]
    P_cot_data.append(fit_phase_shift.results["P_2_prime"].median)
    P_cot_data.append(fit_phase_shift.results["P_2_prime"].em)
    P_cot_data.append(fit_phase_shift.results["P_2_prime"].ep)
    P_cot_data.append(fit_phase_shift.results["P_cot_PS_prime"].median)
    P_cot_data.append(fit_phase_shift.results["P_cot_PS_prime"].em)
    P_cot_data.append(fit_phase_shift.results["P_cot_PS_prime"].ep)
    if "a"+str(order) in fit_phase_shift.result_names:
        for i in range(1+order//2):
            # print(coeff_strs[i]+str(order))
            coeffs.append(np.transpose(fit_phase_shift.results[coeff_strs[i]+str(order)].sample)[0])
    else:
        return None

    mass_Goldstone = fit_phase_shift.results["mass_Goldstone"].median[0]
    return P_cot_data, coeffs, mass_Goldstone, beta, m

def get_phase_shift_mean_min_max_func(coeffs, order = 2, num_inter = 1000):              # returns functions as P_cot_PS (P) (not P squared!)
    coeffs_trans = np.transpose(coeffs)
    P_arr = np.linspace(0,10,num_inter)
    data_mean = np.zeros(num_inter)
    data_p = np.zeros(num_inter)
    data_m = np.zeros(num_inter)

    for i in range(len(P_arr)):
        values_tmp = []
        for j in range(len(coeffs_trans)):
            values_tmp.append(P_cot_delta(P_arr[i], coeffs_trans[j]))
        values_tmp.sort()
        num = len(values_tmp)
        percentage_std = 0.682689
        low_ind = math.ceil(num*(1-percentage_std)/2)
        high_ind = math.floor(num*(1+percentage_std)/2)
        data_mean[i] = values_tmp[num//2]
        data_p[i] = values_tmp[high_ind]
        data_m[i] = values_tmp[low_ind]
    func_mean = inter(x=P_arr, y=data_mean, kind="linear", fill_value="extrapolate")
    func_p = inter(x=P_arr, y=data_p, kind="linear", fill_value="extrapolate")
    func_m = inter(x=P_arr, y=data_m, kind="linear", fill_value="extrapolate")
    return func_mean, func_p, func_m

# def plot_one_P_cot_fit(file, order, color, show_safe = False):
#     data = get_data_from_file(file=file, order=order)
#     if data == None:
#         return None
#     else:
#         P_cot_data, coeffs, mass_Goldstone, beta, m = data
#     plt.errorbar(x=P_cot_data[0], xerr = (P_cot_data[1],P_cot_data[2]), y=P_cot_data[3], yerr=(P_cot_data[4],P_cot_data[5]), color = color, ls = "", capsize=5, markersize=10, label = "b%1.3f, m%1.3f"%(beta,m))
#     func_mean, func_p, func_m = get_phase_shift_mean_min_max_func(coeffs=coeffs, order=order)
#     xarr = np.linspace(0,1)
#     yarr = np.zeros(len(xarr))
#     yarr_p = np.zeros(len(xarr))
#     yarr_m = np.zeros(len(xarr))
#     for i in range(len(xarr)):
#         yarr[i] = func_mean(np.sqrt(xarr[i]))
#         yarr_p[i] = func_p(np.sqrt(xarr[i]))
#         yarr_m[i] = func_m(np.sqrt(xarr[i]))
#     plt.plot(xarr, yarr, color = color)
#     plt.fill_between(xarr, yarr_m, yarr_p, color = color, alpha = 0.1)

#     # plt.xlim((0,1.2))
#     plt.ylim((-4,1))
#     plt.xlabel("$P^{\prime 2}$")
#     plt.ylabel("$P^{\prime}\cot(\delta)$")
#     plt.legend()
#     if show_safe:
#         plt.savefig("plots/P_cot_PS_b%1.3f_m%1.3f_o%i.pdf"%(beta, m, order))
#         plt.title("beta = %1.3f, m1/2 = %1.3f, order = %i"%(beta, m, order))
#         plt.grid()
#         plt.show()
#         plt.clf()

def calc_one_P_cot_fit(file, order):
    data = get_data_from_file(file=file, order=order)
    if data == None:
        return None
    else:
        P_cot_data, coeffs, mass_Goldstone, beta, m = data
    func_mean, func_p, func_m = get_phase_shift_mean_min_max_func(coeffs=coeffs, order=order)
    xarr = np.linspace(0,1)
    yarr = np.zeros(len(xarr))
    yarr_p = np.zeros(len(xarr))
    yarr_m = np.zeros(len(xarr))
    for i in range(len(xarr)):
        yarr[i] = func_mean(np.sqrt(xarr[i]))
        yarr_p[i] = func_p(np.sqrt(xarr[i]))
        yarr_m[i] = func_m(np.sqrt(xarr[i]))
    plot_data = {}
    plot_data["P_cot_data"] = P_cot_data
    plot_data["beta"] = beta
    plot_data["m"] = m
    plot_data["mass_Goldstone"] = mass_Goldstone
    plot_data["P_2"] = xarr
    plot_data["P_cot_PS"] = yarr
    plot_data["P_cot_PS_p"] = yarr_p
    plot_data["P_cot_PS_m"] = yarr_m
    create_plot_file(file+"_P_cot_PS_data_o%i"%order, plot_data)

def plot_one_P_cot_fit(file, order, color, show_safe = False):
    data = read_plot_file(file+"_P_cot_PS_data_o%i"%order)
    if data == None:
        return None
    P_cot_data = data["P_cot_data"]

    beta = data["beta"]
    m = data["m"]
    xarr = data["P_2"]
    yarr = data["P_cot_PS"]
    yarr_p = data["P_cot_PS_p"]
    yarr_m = data["P_cot_PS_m"]
    
    plt.errorbar(x=P_cot_data[0], xerr = (P_cot_data[1],P_cot_data[2]), y=P_cot_data[3], yerr=(P_cot_data[4],P_cot_data[5]), color = color, ls = "", capsize=5, markersize=10, label = "b%1.3f, m%1.3f"%(beta,m))
    plt.plot(xarr, yarr, color = color)
    plt.fill_between(xarr, yarr_m, yarr_p, color = color, alpha = 0.1)

    # plt.xlim((0,1.2))
    plt.ylim((-4,1))
    plt.xlabel("$P^{\prime 2}$")
    plt.ylabel("$P^{\prime}\cot(\delta)$")
    plt.legend()
    if show_safe:
        plt.savefig("plots/P_cot_PS_b%1.3f_m%1.3f_o%i.pdf"%(beta, m, order))
        plt.title("beta = %1.3f, m1/2 = %1.3f, order = %i"%(beta, m, order))
        plt.grid()
        plt.show()
        plt.clf()

def plot_all_P_cot(filenames, orders):
    for i, file in zip(range(len(filenames)), filenames):
        plot_one_P_cot_fit(file, order = orders[i], color = color_arr[i], show_safe=False)
    plt.grid()
    plt.savefig("plots/P_cot_PS.pdf")
    plt.title("universal treshold expansion")
    plt.show()

def E_dispersion_relation(P):       # care for correct use of lattice disp rel, takes P = P/mass_Goldstone, return E = E/mass_Goldstone
    return 2*np.sqrt(1+P**2)

def sigma_of_P(P, P_cot):                                             
    E = E_dispersion_relation(P)
    return 16*np.pi/(E**2*(1+(P_cot/P)**2))

def P_of_s(s):         
    return np.sqrt((s/2)**2-1)

def s_of_P(P):
    return E_dispersion_relation(P)**2

def P_of_v(v, mass_Goldstone):                              # all real unitless
    return v*mass_Goldstone/np.sqrt(1-v*v)

def v_of_P(P, mass_Goldstone):
    return 1/np.sqrt(1+(mass_Goldstone/P)**2)

def calc_one_sigma(file, order):
    data = get_data_from_file(file=file, order=order)
    if data == None:
        return None
    else:
        P_cot_data, coeffs, mass_Goldstone, beta, m = data
    func_mean, func_p, func_m = get_phase_shift_mean_min_max_func(coeffs=coeffs, order=order)

    num = 1000
    # P_max = P_of_s((4*mass_Goldstone)**2, mass_Goldstone)
    # Parr = np.linspace(0,P_max,num)
    Parr = np.linspace(0,3,num)                                                         # P here = P/mass_Goldstone
    sarr = np.linspace(0,10,num)
    varr = np.linspace(0,1,num)
    for i in range(len(Parr)):
        sarr[i] = s_of_P(Parr[i])
        varr[i] = v_of_P(Parr[i], mass_Goldstone=mass_Goldstone)
        # Parr[i] = P_of_v(varr[i], mass_Goldstone)
        # sarr[i] = s_of_P(P_of_v(varr[i], mass_Goldstone), mass_Goldstone)

    sigma_arr = np.zeros(num)
    sigma_arr_p = np.zeros(num)
    sigma_arr_m = np.zeros(num)

    for i in range(len(Parr)):
        sigma_arr[i] = sigma_of_P(Parr[i], func_mean(Parr[i]))
        sigma_arr_p[i] = sigma_of_P(Parr[i], func_p(Parr[i]))
        sigma_arr_m[i] = sigma_of_P(Parr[i], func_m(Parr[i]))
    
    plot_data = {}
    plot_data["P_arr"] = Parr
    plot_data["sarr"] = sarr
    plot_data["varr"] = varr
    plot_data["sigma_arr"] = sigma_arr
    plot_data["sigma_arr_p"] = sigma_arr_p
    plot_data["sigma_arr_m"] = sigma_arr_m
    plot_data["beta"] = beta
    plot_data["m"] = m

    create_plot_file(file+"_sigma_data_o%i"%order, plot_data)
    

    # if xaxis == "P":
    #     plt.plot(Parr, sigma_arr, color = color, label = "b%1.3f m%1.3f"%(beta, m))
    #     plt.fill_between(Parr, sigma_arr_p, sigma_arr_m, color = color, alpha = 0.1)
    # elif xaxis == "v":
    #     plt.plot(varr, sigma_arr, color = color, label = "b%1.3f m%1.3f"%(beta, m))
    #     plt.fill_between(varr, sigma_arr_p, sigma_arr_m, color = color, alpha = 0.1)
    # elif xaxis == "s":
    #     plt.plot(sarr, sigma_arr, color = color, label = "b%1.3f m%1.3f"%(beta, m))
    #     plt.fill_between(sarr, sigma_arr_p, sigma_arr_m, color = color, alpha = 0.1)
    # plt.show()

def plot_one_sigma(file, order, color, show_safe = False, xaxis = "P"):
    # data = get_data_from_file(file=file, order=order)
    # if data == None:
    #     return None
    # else:
    #     P_cot_data, coeffs, mass_Goldstone, beta, m = data
    # func_mean, func_p, func_m = get_phase_shift_mean_min_max_func(coeffs=coeffs, order=order)

    # num = 1000
    # # P_max = P_of_s((4*mass_Goldstone)**2, mass_Goldstone)
    # # Parr = np.linspace(0,P_max,num)
    # Parr = np.linspace(0,3,num)                                                         # P here = P/mass_Goldstone
    # sarr = np.linspace(0,10,num)
    # varr = np.linspace(0,1,num)
    # for i in range(len(Parr)):
    #     sarr[i] = s_of_P(Parr[i])
    #     varr[i] = v_of_P(Parr[i], mass_Goldstone=mass_Goldstone)
    #     # Parr[i] = P_of_v(varr[i], mass_Goldstone)
    #     # sarr[i] = s_of_P(P_of_v(varr[i], mass_Goldstone), mass_Goldstone)

    # sigma_arr = np.zeros(num)
    # sigma_arr_p = np.zeros(num)
    # sigma_arr_m = np.zeros(num)

    plot_data = read_plot_file(file+"_sigma_data_o%i"%order)
    Parr = plot_data["P_arr"]
    sarr = plot_data["sarr"]
    varr = plot_data["varr"]
    sigma_arr = plot_data["sigma_arr"]
    sigma_arr_p = plot_data["sigma_arr_p"]
    sigma_arr_m = plot_data["sigma_arr_m"]
    beta = plot_data["beta"]
    m = plot_data["m"]


    # for i in range(len(Parr)):
    #     sigma_arr[i] = sigma_of_P(Parr[i], func_mean(Parr[i]))
    #     sigma_arr_p[i] = sigma_of_P(Parr[i], func_p(Parr[i]))
    #     sigma_arr_m[i] = sigma_of_P(Parr[i], func_m(Parr[i]))

    if xaxis == "P":
        plt.plot(Parr, sigma_arr, color = color, label = "b%1.3f m%1.3f"%(beta, m))
        plt.fill_between(Parr, sigma_arr_p, sigma_arr_m, color = color, alpha = 0.1)
    elif xaxis == "v":
        plt.plot(varr, sigma_arr, color = color, label = "b%1.3f m%1.3f"%(beta, m))
        plt.fill_between(varr, sigma_arr_p, sigma_arr_m, color = color, alpha = 0.1)
    elif xaxis == "s":
        plt.plot(sarr, sigma_arr, color = color, label = "b%1.3f m%1.3f"%(beta, m))
        plt.fill_between(sarr, sigma_arr_p, sigma_arr_m, color = color, alpha = 0.1)
    # plt.show()

def plot_all_sigma(filenames, orders, xaxis = "P"):
    for i, file in zip(range(len(filenames)), filenames):
        plot_one_sigma(file, color = color_arr[i], xaxis=xaxis, order = orders[i])
    plt.grid()
    plt.ylabel("$\sigma$")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    if xaxis == "P":
        plt.xlabel("P")
        plt.savefig("plots/sigma_of_P.pdf", bbox_inches = "tight")
    if xaxis == "v":
        plt.xlabel("v")
        plt.savefig("plots/sigma_of_v.pdf", bbox_inches = "tight")
    if xaxis == "s":
        plt.xlabel("s")
        plt.savefig("plots/sigma_of_s.pdf", bbox_inches = "tight")
    plt.show()

def MB_dist(v, vmean):
    return 32*v**2*np.exp(-4*v**2/(np.pi*vmean**2))/(np.pi**2*vmean**3)

def intfunc(v, vmean, mass, P_cot_func):
    P = P_of_v(v, mass)
    P_cot = P_cot_func(P)
    return v*sigma_of_P(P, P_cot)*MB_dist(v, vmean)

def plot_integrand(file, color, order = 2):
    data = get_data_from_file(file=file, order=order)
    if data == None:
        return None
    else:
        P_cot_data, coeffs, mass_Goldstone, beta, m = data
    func_mean, func_p, func_m = get_phase_shift_mean_min_max_func(coeffs=coeffs, order=order)
    
    num1 = 1000
    num2 = 20
    # varr = np.linspace(0,1,num1)
    # vmeanarr = np.linspace(0.1,0.9,num2)
    varr = np.logspace(-5,0,num1)
    # vmeanarr = np.logspace(-4,-0.7,num2)
    vmeanarr = np.logspace(-3,-2,num2)
    plt.xscale("log")
    plt.yscale("log")
    sigma_arr =  np.zeros(num1)
    MB_dist_arr = np.zeros((num2,num1))
    intfunc_arr = np.zeros((num2,num1))
    # print(varr)
    print(vmeanarr)
    varr_units = varr*lightspeed

    v_min = 0.1
    v_max = 0.9

    for i in range(num2):
        for j in range(num1):
            MB_dist_arr[i][j] = MB_dist(varr[j], vmeanarr[i])
            intfunc_arr[i][j] = intfunc(varr[j], vmeanarr[i], mass_Goldstone, P_cot_func=func_mean)
        if i == 0:
            plt.plot(varr_units, MB_dist_arr[i], label = "Maxwell_Boltzmann", ls = "--", color = color_arr[i])
            # plt.fill_between(varr_units, intfunc_arr[i], label = "int_func", alpha = 0.3, color = color_arr[i])
        else:
            plt.plot(varr_units, MB_dist_arr[i], ls = "--", color = color_arr[i])
            # plt.fill_between(varr_units, intfunc_arr[i], alpha = 0.3, color = color_arr[i])
            # plt.plot(varr_units, MB_dist_arr[i], label = "MB (v=%1.2f)"%vmeanarr[i], ls = "--", color = color_arr[i])
        plt.fill_between(varr_units, intfunc_arr[i], label = "int_func (v=%e)"%(vmeanarr[i]*lightspeed), color = color_arr[i], alpha = 0.3)
    for i in range(num1):
        P = P_of_v(varr[i], mass_Goldstone)
        P_cot = func_mean(P)
        sigma_arr[i] = sigma_of_P(P, P_cot)
    plt.plot(varr_units, sigma_arr, label = "cross-section", color = "black", linewidth = 2)
    plt.ylim([min(sigma_arr)/10, max(MB_dist_arr[0])*10])

    plt.xlabel("v in km/s")
    plt.title("b%1.3f m%1.3f"%(beta, m))
    plt.grid()
    plt.legend()
    plt.savefig("plots/integrands_b_%1.3f m_%1.3f"%(beta, m)+".pdf", bbox_inches = "tight")
    plt.show()
    plt.clf()

def plot_all_integrands(filenames, orders):
    for i, file in zip(range(len(filenames)), filenames):
        plot_integrand(file, color = color_arr[i], order = orders[i])

def sigma_v_of_vmean(vmean, P_cot_func, mass):
    return scpint.quad(func=intfunc, a=0.0001, b=0.99999, args=(vmean,mass,P_cot_func))[0]    

# def plot_one_sigma_v_vs_v(file, color, order = 2, massscale = 0):
#     data = get_data_from_file(file=file, order=order)
#     if data == None:
#         return None
#     else:
#         P_cot_data, coeffs, mass_Goldstone, beta, m = data
#     func_mean, func_p, func_m = get_phase_shift_mean_min_max_func(coeffs=coeffs, order=order)
#     num = 100
#     # vmeanarr = np.linspace(0, 1, num)
#     vmeanarr = np.logspace(-5, 1, num)
#     sigma_v_arr = np.zeros(num)
#     sigma_v_arr_p = np.zeros(num)
#     sigma_v_arr_m = np.zeros(num)

#     for i in range(len(vmeanarr)):
#         sigma_v_arr[i] = sigma_v_of_vmean(vmeanarr[i], func_mean, mass_Goldstone)/mass_Goldstone
#         sigma_v_arr_p[i] = sigma_v_of_vmean(vmeanarr[i], func_p, mass_Goldstone)/mass_Goldstone
#         sigma_v_arr_m[i] = sigma_v_of_vmean(vmeanarr[i], func_m, mass_Goldstone)/mass_Goldstone
#         # tmp = sigma_v_of_vmean(vmeanarr[i], coeffs, mass_Goldstone)
#         # if tmp < 1e3:
#         #     sigma_v_arr[i] = tmp

#     if massscale != 0:
#         for i in range(len(vmeanarr)):
#             vmeanarr[i] = vmeanarr[i] * lightspeed
#             factor = 6.54692055938448*10**10                                           # 1/MeV^3 in cm^2/g*km/s
#             sigma_v_arr[i] = sigma_v_arr[i]*factor/massscale**3
#             sigma_v_arr_p[i] = sigma_v_arr_p[i]*factor/massscale**3
#             sigma_v_arr_m[i] = sigma_v_arr_m[i]*factor/massscale**3


#     plt.plot(vmeanarr, sigma_v_arr, color = color, label = "b%1.3f m%1.3f"%(beta, m))
#     plt.fill_between(vmeanarr, sigma_v_arr_m, sigma_v_arr_p, color = color, alpha = 0.1)
#     # plt.show()

def calc_one_sigma_v_vs_v(file, order = 2, num = 1000):
    data = get_data_from_file(file=file, order=order)
    if data == None:
        return None
    else:
        P_cot_data, coeffs, mass_Goldstone, beta, m = data
    func_mean, func_p, func_m = get_phase_shift_mean_min_max_func(coeffs=coeffs, order=order)
    # v_mean_arr = np.linspace(0, 1, num)
    v_mean_arr = np.logspace(-4, -2, num)
    sigma_v_arr = np.zeros(num)
    sigma_v_arr_p = np.zeros(num)
    sigma_v_arr_m = np.zeros(num)

    for i in range(len(v_mean_arr)):
        sigma_v_arr[i] = sigma_v_of_vmean(v_mean_arr[i], func_mean, mass_Goldstone)/mass_Goldstone
        sigma_v_arr_p[i] = sigma_v_of_vmean(v_mean_arr[i], func_p, mass_Goldstone)/mass_Goldstone
        sigma_v_arr_m[i] = sigma_v_of_vmean(v_mean_arr[i], func_m, mass_Goldstone)/mass_Goldstone

    plot_data = {}
    plot_data["v_mean_arr"] = v_mean_arr
    plot_data["sigma_v_arr"] = sigma_v_arr
    plot_data["sigma_v_arr_p"] = sigma_v_arr_p
    plot_data["sigma_v_arr_m"] = sigma_v_arr_m
    plot_data["beta"] = beta
    plot_data["m"] = m
    plot_data["mass_Goldstone"] = mass_Goldstone
    
    create_plot_file(file+"_sigma_v_data_o%i"%order, plot_data)

def plot_one_sigma_v_vs_v(file, color, order = 2, massscale = 0):
    data = read_plot_file(file+"_sigma_v_data_o%i"%order)
    v_mean_arr, sigma_v_arr, sigma_v_arr_m, sigma_v_arr_p = [data["v_mean_arr"],data["sigma_v_arr"],data["sigma_v_arr_m"],data["sigma_v_arr_p"]]
    if massscale != 0:
        for i in range(len(v_mean_arr)):
            v_mean_arr[i] = v_mean_arr[i] * lightspeed
            factor = 6.54692055938448*10**10                                           # 1/MeV^3 in cm^2/g*km/s
            sigma_v_arr[i] = sigma_v_arr[i]*factor/massscale**3
            sigma_v_arr_p[i] = sigma_v_arr_p[i]*factor/massscale**3
            sigma_v_arr_m[i] = sigma_v_arr_m[i]*factor/massscale**3
    plt.plot(v_mean_arr, sigma_v_arr, color = color, label = "b%1.3f m%1.3f"%(data["beta"], data["m"]))
    plt.fill_between(v_mean_arr, sigma_v_arr_m, sigma_v_arr_p, color = color, alpha = 0.1)


def plot_all_sigma_v(filenames, orders, massscale = 0):
    for i, file in zip(range(len(filenames)), filenames):
        plot_one_sigma_v_vs_v(file, color = color_arr[i], order = orders[i], massscale=massscale)
    plt.grid()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xscale("log")
    plt.yscale("log")
    if massscale == 0:
        plt.xlabel("$<v>$")
        plt.ylabel("$<\sigma v>$/m")
        plt.savefig("plots/sigma_v_vmean_log_mDM_unitless_zoom.pdf", bbox_inches = "tight")
        plt.show()
    else: 
        plt.title("$m_{DM}$ = %i MeV"%int(massscale))
        HALOdata = np.transpose(np.genfromtxt("input/DMhalodata.dat"))
        plt.scatter(HALOdata[0], HALOdata[1], color = "grey", label = "Halo data")
        plt.xlabel("$<v>$ in km/s")
        plt.ylabel("$<\sigma v>$/m in cm^2/g km/s")
        # plt.legend()
        plt.xlim([2e1,4e3])
        plt.ylim([1e0,1e4])
        plt.savefig("plots/sigma_v_vmean_log_mDM_%i_zoom.pdf"%int(massscale), bbox_inches = "tight")
        plt.xlim([1e1,300000])
        plt.ylim([1e0,1e12])
        plt.savefig("plots/sigma_v_vmean_log_mDM_%i.pdf"%int(massscale), bbox_inches = "tight")
        plt.show()

def main_calc():
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders     = [2,0,0,2,2,4,2,2,4,0,0,0]
    all_order = [0,2,4,6]
    # orders_alt = [2,0,0,0,2,2,2,2,2,0,0,0]    
    # "phase_shift_fit_P_cot_PS_SU(3)_beta5.400_m1-0.890_lim_0.9"

    # for i in range(len(filenames)):
    #     for j in range(len(all_order)):
    #         # calc_one_P_cot_fit(filenames[i], order = all_order[j])
    #         calc_one_sigma(filenames[i], order = all_order[j])
    #         # calc_one_integrand(filenames[i], order = all_order[j])
    for i, file in zip(range(len(filenames)), filenames):
        calc_one_sigma_v_vs_v(file, order = orders[i])
def main_plot():
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders     = [2,0,0,2,2,4,2,2,4,0,0,0]
    all_order = [0,2,4,6]
    # orders_alt = [2,0,0,0,2,2,2,2,2,0,0,0]    
    # "phase_shift_fit_P_cot_PS_SU(3)_beta5.400_m1-0.890_lim_0.9"
    # for i in range(len(filenames)):
    #     for j in range(len(all_order)):
    #         plot_one_P_cot_fit(filenames[i], order = all_order[j], color = color_arr[i], show_safe=True)
    # for i, file in zip(range(len(filenames)), filenames):
    #     plot_one_P_cot_fit(file, order = orders[i], color = color_arr[i], show_safe=True)
    # plot_all_P_cot(filenames, orders)
    # plot_all_sigma(filenames = filenames, orders = orders, xaxis="P")
    # plot_all_sigma(filenames = filenames, orders = orders, xaxis="v")
    # plot_all_sigma(filenames = filenames, orders = orders, xaxis="s")
    # plot_all_integrands(filenames = filenames, orders = orders)
    plot_all_sigma_v(filenames = filenames, orders = orders, massscale=0)
    # plot_all_sigma_v(filenames = filenames, orders = orders, massscale=1)
    # plot_all_sigma_v(filenames = filenames, orders = orders, massscale=10)
    # plot_all_sigma_v(filenames = filenames, orders = orders, massscale=100)

if __name__ == "__main__":
    # main_calc()
    main_plot()