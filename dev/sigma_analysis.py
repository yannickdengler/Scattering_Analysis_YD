import matplotlib.pyplot as plt
import error_classes as errcl 
import numpy as np
import scipy.integrate as scpint
from scipy.interpolate import interp1d as inter
import math
import h5py
import os
import universal_threshold_expansion as UTE
from scipy.special import kn                                # modified Bessel function
from scipy.optimize import bisect
# import mpl_toolkits.mplot3d.axes3d as axes3d

lightspeed = 299792.458
color_arr = ["blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick","blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick","blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick",]

def create_hdf5_file(name, data):
    with h5py.File("/home/dengler_yannick/Documents/Scattering_Analysis_YD/output/plot_files/"+name+".hdf5","w") as f:
        for key, val in data.items():
            f.create_dataset(key, data = val)

def read_hdf5_file(name):
    data = {}
    filename = "/home/dengler_yannick/Documents/Scattering_Analysis_YD/output/plot_files/"+name+".hdf5"
    if os.path.exists(filename):
        with h5py.File(filename,"r") as f:
            for key in f.keys():
                # print(type(f[key]))
                data[key] = f[key][()]
        return data

def v_of_P(P, mass_Goldstone):
    return P/np.sqrt(mass_Goldstone**2+P**2)

def P_of_v(v, mass_Goldstone):
    return v*mass_Goldstone/np.sqrt(1-v*v)

def E_dispersion_relation(P):                           # really 2 Pion energy?
    return 2*np.sqrt(1+P**2)

def s_of_P(P):
    return E_dispersion_relation(P)**2

def M2_of_P(P, P_cot_PS):
    return (16*np.pi)**2/(1+P_cot_PS**2/P**2)

def sigma_of_P(P, P_cot_PS):
    E = E_dispersion_relation(P)
    # print(16*np.pi/(E**2*(1+P_cot_PS**2/P**2)), M2_of_P(P, P_cot_PS)/(16*np.pi*E**2), "should be same")
    return 16*np.pi/(E**2*(1+(P_cot_PS/P)**2))
    # return M2_of_P(P, P_cot_PS)/(16*np.pi*E**2)

def MB_dist(v, vmean):
    return 32*v**2*np.exp(-4*v**2/(np.pi*vmean**2))/(np.pi**2*vmean**3)

def zero_func_vmean_MJ(theta, vmean):
    return 2*theta*(theta+1)*np.exp(-1/theta)/kn(2,1/theta) - vmean

def zero_func_vmean_MJ_exp(theta, vmean):
    sq_t = np.sqrt(theta)
    return sq_t*np.sqrt(8/np.pi)-sq_t**3*7/np.sqrt(8*np.pi)+sq_t**5*105/(32*np.sqrt(2*np.pi))-sq_t**7*525/(256*np.sqrt(2*np.pi))-sq_t**9*9765/(8192*np.sqrt(2*np.pi))-vmean

def gamma_of_v(v):
    return 1/np.sqrt(1-v**2)

def v_of_gamma(gamma):
    return np.sqrt(1-1/gamma**2)

def derivative_gamma_v(v):
    return v/((1-v**2)**(3/2))

def theta_MJ_of_vmean(vmean):
    if vmean < 1e-10:
        return 1e-8
    elif vmean > 0.999999999999997:
        return 1e7
    elif vmean < 0.478:
        return bisect(zero_func_vmean_MJ_exp, 0, 0.11, args=[vmean,])
    else:
        return bisect(zero_func_vmean_MJ, 0.09, 3e7, args=[vmean,])

def MJ_dist(v, vmean):
    return MJ_dist_gamma(gamma_of_v(v), theta_MJ_of_vmean(vmean))

def MJ_dist_gamma(gamma, theta):
    return gamma*np.sqrt(gamma**2-1)*np.exp(-gamma/theta)/(theta*kn(2,1/theta))

def intfunc_MB(v, vmean, mass_Goldstone, P_cot_PS_func):
    P = P_of_v(v, mass_Goldstone)
    P_cot_PS = P_cot_PS_func(P**2)
    return v*sigma_of_P(P, P_cot_PS)*MB_dist(v, vmean)

def sigma_v_of_vmean_MB(vmean, mass_Goldstone, P_cot_PS_func):
    return scpint.quad(intfunc_MB, 0, 1, args=(vmean, mass_Goldstone, P_cot_PS_func))[0]

def intfunc_MJ(gamma, theta, mass_Goldstone, P_cot_PS_func):
    v = v_of_gamma(gamma)
    P = P_of_v(v, mass_Goldstone)
    P_cot_PS = P_cot_PS_func(P**2)
    return v*sigma_of_P(P, P_cot_PS)*MJ_dist_gamma(gamma, theta)

def sigma_v_of_vmean_MJ(vmean, mass_Goldstone, P_cot_PS_func):
    theta = theta_MJ_of_vmean(vmean)
    return scpint.quad(intfunc_MJ, 1, 1e3, args=(theta, mass_Goldstone, P_cot_PS_func))[0]        # check upper limit for gamma (1000 => v=0.9999999)

def calc_distribtions(num_v = 1000, vmean_arr = [0.1,0.5,0.9,0.99]):
    results = {}
    v_arr = np.linspace(0,1,num_v)
    results["v_arr_dist"] = v_arr
    results["gamma_arr_dist"] = gamma_of_v(v_arr)
    for vmean in vmean_arr:
        results["MB_dist_arr_vmean_%f"%vmean] = MB_dist(v_arr, vmean)
        results["MJ_dist_arr_vmean_%f"%vmean] = MJ_dist(v_arr, vmean)
    return results

def calc_M2_sigma_intfunc(file, order, num_P = 1000, vmean_arr = [0.1,0.5,0.9,0.99], update = True):
    results = {}
    if update:
        for key, val in read_hdf5_file(file+"_sigma_analysis").items():
            results[key] = val
    P_cot_data, coeffs_sample, mass_Goldstone, beta, m, N_Ls = UTE.get_data_from_file(file=file, order=order)
    func_mean, func_m, func_p = UTE.get_functions_from_file(file, order)
    # P_2_arr = np.linspace(0, 2*max(P_cot_data[0]), num_P)
    P_2_arr = np.linspace(0, 9, num_P)
    P_arr = np.sqrt(P_2_arr)
    P_cot_PS_mean_arr = func_mean(P_2_arr)
    P_cot_PS_m_arr = func_m(P_2_arr)
    P_cot_PS_p_arr = func_p(P_2_arr)
    v_arr = v_of_P(P_arr, mass_Goldstone)
    gamma_arr = gamma_of_v(v_arr)
    results["P_2_arr"] = P_2_arr
    results["P_arr"] = P_arr
    results["v_arr"] = v_arr
    results["gamma_arr"] = gamma_arr
    results["s_arr"] = s_of_P(P_arr)
    results["P_cot_PS_mean_arr"] = P_cot_PS_mean_arr
    results["P_cot_PS_m_arr"] = P_cot_PS_m_arr
    results["P_cot_PS_p_arr"] = P_cot_PS_p_arr
    M_2_mean_arr, M_2_m_arr, M_2_p_arr, sigma_mean_arr, sigma_m_arr, sigma_p_arr = np.zeros((6,num_P))
    int_func_MB_mean, int_func_MB_m, int_func_MB_p, int_func_MJ_mean, int_func_MJ_m, int_func_MJ_p = np.zeros((6,len(vmean_arr),num_P))
    for i in range(num_P):
        M_2_mean_arr[i] = M2_of_P(P_arr[i], P_cot_PS_mean_arr[i])
        M_2_m_arr[i] = M2_of_P(P_arr[i], P_cot_PS_m_arr[i])
        M_2_p_arr[i] = M2_of_P(P_arr[i], P_cot_PS_p_arr[i])
        sigma_mean_arr[i] = sigma_of_P(P_arr[i], P_cot_PS_mean_arr[i])
        sigma_m_arr[i] = sigma_of_P(P_arr[i], P_cot_PS_m_arr[i])
        sigma_p_arr[i] = sigma_of_P(P_arr[i], P_cot_PS_p_arr[i])
    results["M_2_mean_arr"] = M_2_mean_arr
    results["M_2_m_arr"] = M_2_m_arr
    results["M_2_p_arr"] = M_2_p_arr
    results["sigma_mean_arr"] = sigma_mean_arr
    results["sigma_m_arr"] = sigma_m_arr
    results["sigma_p_arr"] = sigma_p_arr
    for vmean in vmean_arr:
        results["int_func_MB_vmean_%f_mean"%vmean] = intfunc_MB(v_arr, vmean, mass_Goldstone, func_mean)
        results["int_func_MB_vmean_%f_m"%vmean] = intfunc_MB(v_arr, vmean, mass_Goldstone, func_m)
        results["int_func_MB_vmean_%f_p"%vmean] = intfunc_MB(v_arr, vmean, mass_Goldstone, func_p)
        results["int_func_MJ_vmean_%f_mean"%vmean] = intfunc_MJ(gamma_arr, theta_MJ_of_vmean(vmean), mass_Goldstone, func_mean)
        results["int_func_MJ_vmean_%f_m"%vmean] = intfunc_MJ(gamma_arr, theta_MJ_of_vmean(vmean), mass_Goldstone, func_m)
        results["int_func_MJ_vmean_%f_p"%vmean] = intfunc_MJ(gamma_arr, theta_MJ_of_vmean(vmean), mass_Goldstone, func_p)
    return results

def calc_sigma_v(file, order, num_v = 1000):
    P_cot_data, coeffs_sample, mass_Goldstone, beta, m, N_Ls = UTE.get_data_from_file(file=file, order=order)
    func_mean, func_m, func_p = UTE.get_functions_from_file(file, order)
    results = {}
    v_mean_arr = np.logspace(-4.5,-0.0001, num_v)
    P_mean_arr = P_of_v(v_mean_arr, mass_Goldstone)
    results["v_mean_arr"] = v_mean_arr
    results["P_mean_arr"] = P_mean_arr
    results["s_mean_arr"] = s_of_P(P_mean_arr)
    sigma_v_MB_mean_arr, sigma_v_MB_m_arr, sigma_v_MB_p_arr, sigma_v_MJ_mean_arr, sigma_v_MJ_m_arr, sigma_v_MJ_p_arr = np.zeros((6,num_v))
    for i in range(num_v):
        sigma_v_MB_mean_arr[i] = sigma_v_of_vmean_MB(v_mean_arr[i], mass_Goldstone, func_mean)
        sigma_v_MB_m_arr[i] = sigma_v_of_vmean_MB(v_mean_arr[i], mass_Goldstone, func_m)
        sigma_v_MB_p_arr[i] = sigma_v_of_vmean_MB(v_mean_arr[i], mass_Goldstone, func_p)
        sigma_v_MJ_mean_arr[i] = sigma_v_of_vmean_MJ(v_mean_arr[i], mass_Goldstone, func_mean)
        sigma_v_MJ_m_arr[i] = sigma_v_of_vmean_MJ(v_mean_arr[i], mass_Goldstone, func_m)
        sigma_v_MJ_p_arr[i] = sigma_v_of_vmean_MJ(v_mean_arr[i], mass_Goldstone, func_p)
    results["sigma_v_MB_mean_arr"] = sigma_v_MB_mean_arr
    results["sigma_v_MB_m_arr"] = sigma_v_MB_m_arr
    results["sigma_v_MB_p_arr"] = sigma_v_MB_p_arr
    results["sigma_v_MJ_mean_arr"] = sigma_v_MJ_mean_arr
    results["sigma_v_MJ_m_arr"] = sigma_v_MJ_m_arr
    results["sigma_v_MJ_p_arr"] = sigma_v_MJ_p_arr
    return results

def calc(file, order):
    results = {}
    # for key, val in calc_distribtions().items():
    #     results[key] = val
    for key, val in calc_M2_sigma_intfunc(file, order).items():
        results[key] = val
    # for key, val in calc_sigma_v(file, order).items():
    #     results[key] = val
    # print(results)
    create_hdf5_file(file+"_sigma_analysis", results)
    print(read_hdf5_file(file+"_sigma_analysis"))


def main_calc():
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders     = [2,0,0,2,2,4,2,2,4,0,0,0]
    for i in range(len(filenames)):
    # for i in range(1):
        calc(filenames[i], orders[i])

def plot_one_sigma_v_vs_v(filename, color, order = 2, massscale = 0, error = True, xaxis = "v_mean", dist = "both"):
    print(filename+"_sigma_analysis")
    data = read_hdf5_file(filename+"_sigma_analysis")
    P_cot_data, coeffs_sample, mass_Goldstone, beta, m, N_Ls = UTE.get_data_from_file(file=filename, order=order)
    x_arr = data["%s_arr"%xaxis]
    y_arr_MB = [data["sigma_v_MB_mean_arr"],data["sigma_v_MB_m_arr"],data["sigma_v_MB_p_arr"]]
    y_arr_MJ = [data["sigma_v_MJ_mean_arr"],data["sigma_v_MJ_m_arr"],data["sigma_v_MJ_p_arr"]]

    if massscale != 0:
        for i in range(len(x_arr)):
            x_arr[i] = x_arr[i] * lightspeed
            for j in range(len(y_arr_MB)):
                y_arr_MB[j][i] = y_arr_MB[j][i]/(massscale**3*4.578*1e-6)*lightspeed
                y_arr_MJ[j][i] = y_arr_MJ[j][i]/(massscale**3*4.578*1e-6)*lightspeed
    if dist == "MB" or dist == "both":
        plt.plot(x_arr, y_arr_MB[0], color = color, label = "b%1.3f m%1.3f_MB"%(beta, m))
    if dist == "MJ" or dist == "both":
        plt.plot(x_arr, y_arr_MJ[0], color = color, label = "b%1.3f m%1.3f_MJ"%(beta, m), ls = "--")
    if error:
        if dist == "MB" or dist == "both":
            plt.fill_between(x_arr, y_arr_MB[1], y_arr_MB[2], color = color, alpha = 0.1)
        elif dist == "MJ" or dist == "both":
            plt.fill_between(x_arr, y_arr_MJ[1], y_arr_MJ[2], color = color, alpha = 0.1)
    # plt.show()

def plot_all_sigma_v(massscale = 0, xaxis = "v_mean"):
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders     = [2,0,0,2,2,4,2,2,4,0,0,0]
    # orders_alt = [2,0,0,0,2,2,2,2,2,0,0,0]    
    # "phase_shift_fit_P_cot_PS_SU(3)_beta5.400_m1-0.890_lim_0.9"
    for i, file in zip(range(len(filenames)), filenames):
        plot_one_sigma_v_vs_v(file, color = color_arr[i], order = orders[i], massscale=massscale)
    # plot_one_sigma_v_vs_v(filenames[0], color = color_arr[0], order = orders[0], massscale=massscale, dist="both")
    plt.grid()
    if massscale == 0:
        plt.xlabel("$<v>$")
        plt.ylabel("$<\sigma v>$/m")
    else: 
        HALOdata = np.transpose(np.genfromtxt("input/DMhalodata.dat"))
        plt.scatter(HALOdata[0], HALOdata[1], color = "grey", label = "Halo data")
        plt.xlabel("$<v>$ in km/s")
        plt.ylabel("$<\sigma v>$/m in cm^2/g km/s")
        v_arr = np.logspace(np.log10(lightspeed)-6,np.log10(lightspeed), 20)
        # print(v_arr)
        for i in range(20):
            cross_section_mass = 10**(i-16)
            # print("hey")
            plt.plot(v_arr,cross_section_mass*v_arr, color = "grey", ls = "--")
    # plt.legend()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # plt.savefig("plots/sigma_v_vmean_mDM_%i.pdf"%int(massscale), bbox_inches = "tight")
    # plt.show()
    plt.title("$m_{DM}$ = %i MeV"%int(massscale))
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("plots/sigma_v_vmean_log_mDM_%i.pdf"%int(massscale), bbox_inches = "tight")
    plt.show()

def plot_M_2(x_axis = "P"):
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders = [2,0,0,2,2,4,2,2,4,0,0,0]
    for i, filename in zip(range(len(filenames)), filenames):
        data = read_hdf5_file(filename+"_sigma_analysis")
        P_cot_data, coeffs_sample, mass_Goldstone, beta, m, N_Ls = UTE.get_data_from_file(file=filename, order=orders[i])
        y_arr_M2 = [data["M_2_mean_arr"],data["M_2_m_arr"],data["M_2_p_arr"]]
        x_arr = data["%s_arr"%x_axis]
        plt.xlabel(x_axis)
        plt.ylabel("$|M|^2$")
        plt.grid()
        plt.plot(x_arr,y_arr_M2[0], color = color_arr[i], label = "b%f, m%f"%(beta, m))
        plt.fill_between(x_arr,y_arr_M2[1],y_arr_M2[2], color = color_arr[i], alpha = 0.1)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    plt.savefig("plots/all_M_2_vs_%s.pdf"%x_axis, bbox_inches = "tight")
    plt.clf()
def plot_sigma(x_axis = "P"):
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders = [2,0,0,2,2,4,2,2,4,0,0,0]
    for i, filename in zip(range(len(filenames)), filenames):
        data = read_hdf5_file(filename+"_sigma_analysis")
        P_cot_data, coeffs_sample, mass_Goldstone, beta, m, N_Ls = UTE.get_data_from_file(file=filename, order=orders[i])
        y_arr_sigma = [data["sigma_mean_arr"],data["sigma_m_arr"],data["sigma_p_arr"]]
        x_arr = data["%s_arr"%x_axis]
        plt.xlabel(x_axis)
        plt.ylabel("$\sigma$")
        plt.grid()
        plt.plot(x_arr,y_arr_sigma[0], color = color_arr[i], label = "b%f, m%f"%(beta, m))
        plt.fill_between(x_arr,y_arr_sigma[1],y_arr_sigma[2], color = color_arr[i], alpha = 0.1)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    plt.savefig("plots/all_sigmas_vs_%s.pdf"%x_axis, bbox_inches = "tight")
    plt.clf()

def plot_intfunc(x_axis = "P"):
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders = [2,0,0,2,2,4,2,2,4,0,0,0]
    for i, filename in zip(range(len(filenames)), filenames):
        data = read_hdf5_file(filename+"_sigma_analysis")
        P_cot_data, coeffs_sample, mass_Goldstone, beta, m, N_Ls = UTE.get_data_from_file(file=filename, order=orders[i])
        x_arr = data["%s_arr"%x_axis]
        plt.xlabel(x_axis)
        plt.ylabel("$\sigma$")
        plt.grid()
        vmean_arr = [0.1,0.5,0.9,0.99]
        for j in range(len(vmean_arr)):
            y_arr_intfunc_MB = [data["int_func_MB_vmean_%f_mean"%vmean_arr[j]],data["int_func_MB_vmean_%f_m"%vmean_arr[j]],data["int_func_MB_vmean_%f_p"%vmean_arr[j]]]
            y_arr_intfunc_MJ = [data["int_func_MJ_vmean_%f_mean"%vmean_arr[j]],data["int_func_MJ_vmean_%f_m"%vmean_arr[j]],data["int_func_MJ_vmean_%f_p"%vmean_arr[j]]]
            plt.title("b%f, m%f"%(beta, m))
            plt.plot(x_arr,y_arr_intfunc_MB[0], color = color_arr[j], label = "vmean = %f"%vmean_arr[j])
            plt.fill_between(x_arr,y_arr_intfunc_MB[1],y_arr_intfunc_MB[2], color = color_arr[j], alpha = 0.1)
            plt.plot(x_arr,y_arr_intfunc_MJ[0], color = color_arr[j], ls="--")
            plt.fill_between(x_arr,y_arr_intfunc_MJ[1],y_arr_intfunc_MJ[2], color = color_arr[j], alpha = 0.1)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig("plots/intfunc_vs_%s_b%s_m%s.pdf"%(x_axis,beta,m), bbox_inches = "tight")
        plt.show()
        plt.clf()

def plot_sigma_v_dist_intfunc(x_axis="v"):
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders = [2,0,0,2,2,4,2,2,4,0,0,0]
    for i, filename in zip(range(len(filenames)), filenames):
        data = read_hdf5_file(filename+"_sigma_analysis")
        P_cot_data, coeffs_sample, mass_Goldstone, beta, m, N_Ls = UTE.get_data_from_file(file=filename, order=orders[i])
        x_arr = data["%s_arr"%x_axis]
        plt.xlabel(x_axis)
        plt.ylabel("$\sigma$")
        plt.grid()
        vmean_arr = [0.1,0.5,0.9,0.99]
        for j in range(len(vmean_arr)):
            plt.plot(x_arr,data["v_arr"], color = "black", label = "v")
            plt.plot(x_arr,data["MB_dist_arr_vmean_%f"%vmean_arr[j]], label = "MB_dist")
            plt.plot(x_arr,data["MJ_dist_arr_vmean_%f"%vmean_arr[j]], label = "MJ_dist", ls = "--")
            plt.plot(x_arr,data["sigma_mean_arr"], label = "$\sigma$")
            plt.plot(x_arr,data["int_func_MB_vmean_%f_mean"%vmean_arr[j]], label = "intfunc_MB")
            plt.plot(x_arr,data["int_func_MJ_vmean_%f_mean"%vmean_arr[j]], label = "intfunc_MJ", ls = "--")
            plt.title("b%f, m%f, vmean%f"%(beta, m, vmean_arr[j]))
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.savefig("plots/intfunc_v_dist_vs_%s_b%s_m%s.pdf"%(x_axis,beta,m), bbox_inches = "tight")
            plt.show()
            plt.clf()

# def a_0_m_pi_vs_L(a_0s, N_Ls, betas, ms):
#     filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
#     orders = [2,0,0,2,2,4,2,2,4,0,0,0]
#     for i, filename in zip(range(len(filenames)), filenames):
#         data = read_hdf5_file(filename+"_sigma_analysis")
#         P_cot_data, coeffs_sample, mass_Goldstone, beta, m, N_Ls = UTE.get_data_from_file(file=filename, order=orders[i])
#         num = len(coeffs_sample[0])
#         percentage_std = 0.682689
#         low_ind = math.ceil(num*(1-percentage_std)/2)
#         high_ind = math.floor(num*(1+percentage_std)/2)
#         coeffs_mean = np.transpose(coeffs_sample)[num//2]
#         coeffs_m = np.transpose(coeffs_sample)[low_ind]
#         coeffs_p = np.transpose(coeffs_sample)[high_ind]
#         print(coeffs_mean, coeffs_m, coeffs_p)
#         a_0 = [1/coeffs_mean[0],(1/coeffs_mean[0]-1/coeffs_p[0]),(1/coeffs_mean[0]-1/coeffs_m[0])]
#         plt.errorbar(x=N_Ls[i], y=a_0s[i][0], yerr=[a_0s[i][1], a_0s[i][2]], color = color_arr[i], ls = "", capsize=5, markersize=10, label = "b%1.3f, m%1.3f"%(betas[i],ms[i]))
#     plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#     plt.savefig("plots/a_0s_vs_L.pdf", bbox_inches = "tight")
#     plt.show()

def a_b_heatmap():
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders = [2,0,0,2,2,4,2,2,4,0,0,0]
    for i, filename in zip(range(len(filenames)), filenames):
        data = read_hdf5_file(filename+"_sigma_analysis")
        P_cot_data, coeffs_sample, mass_Goldstone, beta, m, N_Ls = UTE.get_data_from_file(file=filename, order=orders[i])
        num = len(coeffs_sample[0])
        percentage_std = 0.682689
        low_ind = math.ceil(num*(1-percentage_std)/2)
        high_ind = math.floor(num*(1+percentage_std)/2)
        coeffs_mean = np.transpose(coeffs_sample)[num//2]
        coeffs_m = np.transpose(coeffs_sample)[low_ind]
        coeffs_p = np.transpose(coeffs_sample)[high_ind]
        print(coeffs_mean, coeffs_m, coeffs_p)
        if len(coeffs_mean) == 1:
            plt.errorbar(x=coeffs_mean, y=0, xerr=[abs(coeffs_mean-coeffs_m),abs(coeffs_mean-coeffs_p)], color = color_arr[i], ls = "", capsize=5, markersize=10, label = "b%1.3f, m%1.3f"%(beta,m))
        else:
            plt.errorbar(x=coeffs_mean[0], y=coeffs_mean[1], xerr=[[abs(coeffs_mean[0]-coeffs_m[0]),],[abs(coeffs_mean[0]-coeffs_p[0]),]], yerr=[[abs(coeffs_mean[1]-coeffs_m[1]),],[abs(coeffs_mean[1]-coeffs_p[1]),]], color = color_arr[i], ls = "", capsize=5, markersize=10, label = "b%1.3f, m%1.3f"%(beta,m))
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.grid()
    plt.xlabel("a")
    plt.ylabel("b")
    plt.savefig("plots/a_b_heatmap.pdf", bbox_inches="tight")
    plt.xlim([-6,1])
    plt.ylim([-2,12])
    plt.savefig("plots/a_b_heatmap_zoom.pdf", bbox_inches="tight")
    plt.show()

def re_b_prime(b_prime, m):
    return 2*b_prime/m
def a_a_prime(a_prime, m):
    return -1/(a_prime*m)

def a_b_heatmap_advanced():
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    mass_arr = np.linspace(0,25, 100)
    a_arr = []
    re_arr = []
    orders = [2,0,0,2,2,4,2,2,4,0,0,0]
    betas = []
    ms = []
    for i, filename in zip(range(len(filenames)), filenames):
        a_arr.append([])
        re_arr.append([])
        data = read_hdf5_file(filename+"_sigma_analysis")
        P_cot_data, coeffs_sample, mass_Goldstone, beta, m, N_Ls = UTE.get_data_from_file(file=filename, order=orders[i])
        betas.append(beta)
        ms.append(m)
        num = len(coeffs_sample[0])
        percentage_std = 0.682689
        low_ind = math.ceil(num*(1-percentage_std)/2)
        high_ind = math.floor(num*(1+percentage_std)/2)
        coeffs_mean = np.transpose(coeffs_sample)[num//2]
        # coeffs_m = np.transpose(coeffs_sample)[low_ind]
        # coeffs_p = np.transpose(coeffs_sample)[high_ind]
        MeVfm = 197.3
        for j in range(len(mass_arr)):
            if len(coeffs_mean) == 1:
                # a_arr[i].append(a_a_prime(coeffs_mean, mass_arr[j])*MeVfm)
                a_arr[i].append(a_a_prime(coeffs_mean, mass_arr[j])*1000*MeVfm)
                re_arr[i].append(0)
            else:
                # a_arr[i].append(a_a_prime(coeffs_mean[0], mass_arr[j])*MeVfm)
                # re_arr[i].append(re_b_prime(coeffs_mean[1],mass_arr[j])*MeVfm)
                a_arr[i].append(a_a_prime(coeffs_mean[0], mass_arr[j])*1000*MeVfm)
                re_arr[i].append(re_b_prime(coeffs_mean[1],mass_arr[j])*1000*MeVfm)
    fig, axs = plt.subplots(nrows=2, ncols=2)
    axs[0,1].axis("off")
    axs[0,0].grid()
    axs[0,1].grid()
    axs[1,0].grid()
    axs[1,1].grid()
    # axs[0,0].set_xlabel("")
    axs[0,0].set_ylabel("$r_e$ [fm]")
    axs[1,0].set_xlabel("a [fm]")
    axs[1,0].set_ylabel("$m_{DM}$ [GeV]")
    axs[1,1].set_xlabel("$r_e$ [fm]")
    # axs[1,1].set_ylabel([0,25])
    for i in range(len(a_arr)):
        axs[0,0].plot(a_arr[i], re_arr[i], color = color_arr[i], label = "b%1.3f, m%1.3f"%(betas[i], ms[i]))
        axs[1,0].plot(a_arr[i], mass_arr, color = color_arr[i])
        axs[1,1].plot(re_arr[i], mass_arr, color = color_arr[i])
    axs[0,0].legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol = 2)

    plt.savefig("plots/a_b_heatmap.pdf", bbox_inches = "tight")
    axs[0,0].set_xscale('log')
    axs[0,0].set_yscale('log')
    axs[1,0].set_xscale('log')
    axs[1,0].set_yscale('log')
    axs[1,1].set_xscale('log')
    axs[1,1].set_yscale('log')
    plt.savefig("plots/a_b_heatmap_log.pdf", bbox_inches = "tight")
    axs[0,0].set_xscale('linear')
    axs[0,0].set_yscale('linear')
    axs[1,0].set_xscale('linear')
    axs[1,0].set_yscale('linear')
    axs[1,1].set_xscale('linear')
    axs[1,1].set_yscale('linear')
    axs[0,0].set_xlim([0,40])
    axs[0,0].set_ylim([-100,100])
    axs[1,0].set_xlim([0,40])
    axs[1,0].set_ylim([0,25])
    axs[1,1].set_xlim([-100,100])
    axs[1,1].set_ylim([0,25])
    plt.savefig("plots/a_b_heatmap_zoom.pdf", bbox_inches = "tight")

    # plt.show()


def main_plot():
    a_b_heatmap_advanced()
    # plot_all_sigma_v(massscale=10)
    # plot_M_2()
    # plot_sigma()
    # plot_M_2(x_axis="s")
    # plot_sigma(x_axis="s")
    # plot_M_2(x_axis="v")
    # plot_sigma(x_axis="v")
    # plot_intfunc()
    # plot_intfunc(x_axis="v")
    # plot_sigma_v_dist_intfunc()



if __name__ == "__main__":
    # main_calc()
    main_plot()
    # extra_plots()

