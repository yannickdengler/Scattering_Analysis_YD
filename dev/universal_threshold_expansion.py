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

def P_cot_delta(P2, coeffs):                                 # takes and return real unitless P = P/mass_Goldstone
    res = 0
    for i in range(len(coeffs)):
        if i == 0:
            res += coeffs[i]
        else:
            res += coeffs[i]*P2**(i)
    return res

def get_data_from_file(file, order=2):
    coeff_strs = ["a", "b", "c", "d", "e"]
    fit_phase_shift = errcl.measurement(file)
    fit_phase_shift.read_from_HDF()
    P_cot_data = []
    coeffs = []
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
            coeffs.append(np.transpose(fit_phase_shift.results[coeff_strs[i]+str(order)].sample)[0])
    else:
        return None
    mass_Goldstone = fit_phase_shift.results["mass_Goldstone"].median[0]
    return P_cot_data, coeffs, mass_Goldstone, beta, m, N_Ls

def perform_fit(coeffs_sample, num_inter = 1000):
    coeffs_trans = np.transpose(coeffs_sample)
    P2_inter = np.linspace(0,2,num_inter)
    P_cot_inter = np.zeros(num_inter)
    P_cot_inter_p = np.zeros(num_inter)
    P_cot_inter_m = np.zeros(num_inter)
    for i in range(len(P2_inter)):
        values_tmp = []
        for j in range(len(coeffs_trans)):
            values_tmp.append(P_cot_delta(P2_inter[i], coeffs_trans[j]))
        values_tmp.sort()
        num = len(values_tmp)
        percentage_std = 0.682689
        low_ind = math.ceil(num*(1-percentage_std)/2)
        high_ind = math.floor(num*(1+percentage_std)/2)
        P_cot_inter[i] = values_tmp[num//2]
        P_cot_inter_p[i] = values_tmp[high_ind]
        P_cot_inter_m[i] = values_tmp[low_ind]
    return P2_inter, P_cot_inter, P_cot_inter_m, P_cot_inter_p

def create_function(P2_inter, P_cot_inter, P_cot_inter_m, P_cot_inter_p):
    func_mean = inter(x=P2_inter, y=P_cot_inter, kind="linear", fill_value="extrapolate")
    func_m = inter(x=P2_inter, y=P_cot_inter_m, kind="linear", fill_value="extrapolate")
    func_p = inter(x=P2_inter, y=P_cot_inter_p, kind="linear", fill_value="extrapolate")
    return func_mean, func_m, func_p

def get_functions_from_file(file, order=2):                         # gives P_cot_PS as function of P^2
    P_cot_data, coeffs_sample, mass_Goldstone, beta, m, N_Ls = get_data_from_file(file=file, order=order)
    P2_inter, P_cot_inter, P_cot_inter_m, P_cot_inter_p = perform_fit(coeffs_sample)
    return create_function(P2_inter, P_cot_inter, P_cot_inter_m, P_cot_inter_p)
    # return create_function(perform_fit(coeffs_sample))

def calc_P_cot_fit(file, order=2):
    data = get_data_from_file(file=file, order=order)
    if data == None:
        return None
    else:
        P_cot_data, coeffs_sample, mass_Goldstone, beta, m, N_Ls = data
    plot_data = {}
    plot_data["P2_prime"] = P_cot_data[0]
    plot_data["P2_prime_m"] = P_cot_data[1]
    plot_data["P2_prime_p"] = P_cot_data[2]
    plot_data["P_cot_PS_prime"] = P_cot_data[3]
    plot_data["P_cot_PS_prime_m"] = P_cot_data[4]
    plot_data["P_cot_PS_prime_p"] = P_cot_data[5]
    plot_data["P2_inter"], plot_data["P_cot_inter"], plot_data["P_cot_inter_m"], plot_data["P_cot_inter_p"] = perform_fit(coeffs_sample)
    plot_data["beta"] = beta
    plot_data["m"] = m
    plot_data["mass_Goldstone"] = mass_Goldstone
    create_plot_file(file+"_P_cot_PS_data_o%i"%order, plot_data)

def get_data_from_functions(functions, num_x = 100):
    P2_arr = np.linspace(0,1.2,num_x)
    return P2_arr, functions[0](P2_arr), functions[1](P2_arr), functions[2](P2_arr)

def plot_P_cot_PS(file, order = 2, color = "grey", xaxis = "P", fit = True, scatter = False, show = True, save = True):
    pd = read_plot_file(file+"_P_cot_PS_data_o%i"%order)               # plot_data
    if pd == None:
        return None
    beta, m, mass_Goldstone = [pd["beta"],pd["m"],pd["mass_Goldstone"]]
    if fit:
        fit_data = get_data_from_functions(create_function(pd["P2_inter"],pd["P_cot_inter"],pd["P_cot_inter_m"],pd["P_cot_inter_p"]))
    if xaxis == "P":
        plt.xlim((0,1))
        plt.ylim((-4,1))
        plt.xlabel("$P^{\prime}$")
        plt.ylabel("$P^{\prime}\cot(\delta)$")
        xaxis_data = np.sqrt(pd["P2_prime"])
        xaxis_inter = np.sqrt(pd["P2_inter"])
        xaxis_fit = np.sqrt(fit_data[0])
    if xaxis == "P2":
        plt.xlim((0,0.4))
        plt.ylim((-4,1))
        plt.xlabel("$P^{\prime 2}$")
        plt.ylabel("$P^{\prime}\cot(\delta)$")
        xaxis_data = pd["P2_prime"]
        xaxis_inter = pd["P2_inter"]
        xaxis_fit = fit_data[0]
    if xaxis == "v":
        plt.xlim((0,300000))
        plt.ylim((-4,1))
        plt.xlabel("v in km/s")
        plt.ylabel("$P^{\prime}\cot(\delta)$")
        xaxis_data = v_of_P(pd["P2_prime"], mass_Goldstone)*lightspeed
        xaxis_inter = v_of_P(pd["P2_inter"], mass_Goldstone)*lightspeed
        xaxis_fit = v_of_P(fit_data[0], mass_Goldstone)*lightspeed
    if xaxis == "s":
        plt.xlim((3.9,s_of_P(0.4)))
        plt.ylim((-4,1))
        plt.xlabel("$s^\prime$")
        plt.ylabel("$P^{\prime}\cot(\delta)$")
        xaxis_data = s_of_P(pd["P2_prime"])
        xaxis_inter = s_of_P(pd["P2_inter"])
        xaxis_fit = s_of_P(fit_data[0])
        plt.axvline(4)

    plt.errorbar(x=xaxis_data,xerr=(pd["P2_prime_m"],pd["P2_prime_p"]),y=pd["P_cot_PS_prime"],yerr=(pd["P_cot_PS_prime_m"],pd["P_cot_PS_prime_p"]), color = color, ls = "", capsize=5, markersize=10, label = "b%1.3f, m%1.3f"%(beta,m))
    if scatter:
        plt.scatter(xaxis_inter, pd["P_cot_inter"], color = "black", marker = "x", s = 1)
        plt.scatter(xaxis_inter, pd["P_cot_inter_m"], color = "black", marker = "v", s = 1)
        plt.scatter(xaxis_inter, pd["P_cot_inter_p"], color = "black", marker = "^", s = 1)
    if fit:
        # fit_data = get_data_from_functions(create_function(pd["P2_inter"],pd["P_cot_inter"],pd["P_cot_inter_m"],pd["P_cot_inter_p"]))
        plt.plot(xaxis_fit,fit_data[1], color = color)
        plt.fill_between(x=xaxis_fit,y1=fit_data[2],y2=fit_data[3], color = color, alpha = 0.1)

    if save:
        plt.title("beta = %1.3f, m1/2 = %1.3f, order = %i"%(beta, m, order))
        plt.grid()
        plt.savefig("plots/P_cot_PS_vs_%s_b%1.3f_m%1.3f_o%i.pdf"%(xaxis,beta,m,order))
    if show:
        plt.show()
    if save or show:
        plt.clf()

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

def main_calc():
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders     = [2,0,0,2,2,4,2,2,4,0,0,0]
    all_order = [0,2,4,6]
    # orders_alt = [2,0,0,0,2,2,2,2,2,0,0,0]    
    # "phase_shift_fit_P_cot_PS_SU(3)_beta5.400_m1-0.890_lim_0.9"

    for i in range(len(filenames)):
        for j in range(len(all_order)):
            calc_P_cot_fit(filenames[i], order = all_order[j])
def main_plot():
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders     = [2,0,0,2,2,4,2,2,4,0,0,0]
    all_order = [0,2,4,6]
    # orders_alt = [2,0,0,0,2,2,2,2,2,0,0,0]    
    # "phase_shift_fit_P_cot_PS_SU(3)_beta5.400_m1-0.890_lim_0.9"
    for i in range(len(filenames)):
        for j in range(len(all_order)):
            plot_P_cot_PS(filenames[i], order = all_order[j], color = color_arr[i], xaxis="P", show=True, save=True)
            # plot_P_cot_PS(filenames[i], order = all_order[j], color = color_arr[i], xaxis="P2", show=True, save=True)
            # plot_P_cot_PS(filenames[i], order = all_order[j], color = color_arr[i], xaxis="s", show=True, save=True)
            # plot_P_cot_PS(filenames[i], order = all_order[j], color = color_arr[i], xaxis="v", show=True, save=True)

if __name__ == "__main__":
    # main_calc()
    main_plot()