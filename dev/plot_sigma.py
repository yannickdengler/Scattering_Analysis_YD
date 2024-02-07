import matplotlib.pyplot as plt
import error_classes as errcl 
import numpy as np
import scipy.integrate as scpint
import math
import h5py

lightspeed = 299792.458
color_arr = ["blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick","blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick","blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick",]

def get_data_from_file(beta, m_1, limit = 0.9, order=2):
    coeff_strs = ["a", "b", "c", "d", "e"]
    P_cot_data = []
    coeffs = []
    coeffs_p = []
    coeffs_m = []
    with h5py.File("/home/dengler_yannick/Documents/Scattering_Analysis_YD/output/phase_shift_fit_results/phase_shift_fit_results_b%1.3f_m%1.3f_lim%1.3f"%(beta,m_1, limit),"r") as f:
        P_2_median, P_2_ep, P_2_em, P_cot_PS_median, P_cot_PS_ep, P_cot_PS_em = f["plot_data"][:]
        P_cot_data.append(P_2_median)
        P_cot_data.append(P_2_ep)
        P_cot_data.append(P_cot_PS_median)
        P_cot_data.append(P_cot_PS_ep)
        coeffs = []
        coeffs_p = []
        coeffs_m = []
        for i in range(1+order//2):
            coeffs.append(float(f[coeff_strs[i]+str(order)+"_mean"][:]))
            coeffs_p.append(float(f[coeff_strs[i]+str(order)+"_p"][:]))
            coeffs_m.append(float(f[coeff_strs[i]+str(order)+"_m"][:]))
        mass_Goldstone = float(f["mass_Goldstone_median"][:])
        beta = float(f["beta"][:])
        m = float(f["m_1"][:])

    return P_cot_data, coeffs, coeffs_p, coeffs_m, mass_Goldstone, beta, m

# def get_data_from_file(file, order=2):
#     coeff_strs = ["a", "b", "c", "d", "e"]
#     fit_phase_shift = errcl.measurement(file)
#     fit_phase_shift.read_from_HDF()
#     P_cot_data = []
#     coeffs = []
#     coeffs_p = []
#     coeffs_m = []

#     beta = fit_phase_shift.infos["beta"]
#     m = fit_phase_shift.infos["m_1"]
#     P_cot_data.append(fit_phase_shift.results["P_2_prime"].median)
#     P_cot_data.append(fit_phase_shift.results["P_2_prime"].e)
#     P_cot_data.append(fit_phase_shift.results["P_cot_PS_prime"].median)
#     P_cot_data.append(fit_phase_shift.results["P_cot_PS_prime"].e)
#     if "a"+str(order) in fit_phase_shift.result_names:
#         for i in range(1+order//2):
#             coeffs.append(fit_phase_shift.results[coeff_strs[i]+str(order)].median[0])
#             # coeffs_p.append(fit_phase_shift.results[coeff_strs[i]+str(order)].median_p[0])
#             # coeffs_m.append(fit_phase_shift.results[coeff_strs[i]+str(order)].median_m[0])
#     else:
#         return None
#     val_coeffs_arr = []
#     val_coeffs_arr.append(np.transpose(fit_phase_shift.results["integral"+str(order)].sample)[0])
#     for i in range(1+order//2):
#         val_coeffs_arr.append(np.transpose(fit_phase_shift.results[coeff_strs[i]+str(order)].sample)[0])
#     for i in range(1,len(val_coeffs_arr)):
#         val_coeffs_arr[i] = [x for _, x in sorted(zip(val_coeffs_arr[0],val_coeffs_arr[i]))]
#     val_coeffs_arr[0] = sorted(val_coeffs_arr[0])
#     num = len(val_coeffs_arr[0])
#     percentage_std = 0.682689
#     low_ind = math.ceil(num*(1-percentage_std)/2)
#     high_ind = math.floor(num*(1+percentage_std)/2)
#     for i in range(1,len(val_coeffs_arr)):
#         coeffs_p.append(val_coeffs_arr[i][high_ind])
#         coeffs_m.append(val_coeffs_arr[i][low_ind])

#     mass_Goldstone = fit_phase_shift.results["mass_Goldstone"].median[0]

#     return P_cot_data, coeffs, coeffs_p, coeffs_m, mass_Goldstone, beta, m

def P_cot_delta(P, coeffs):                                 # takes and return real unitless P = P/mass_Goldstone
    res = 0
    for i in range(len(coeffs)):
        if i == 0:
            res += coeffs[i]
        else:
            res += coeffs[i]*P**(2*i)
    return res

def P_cot_delta_min_max(P, coeffs_min, coeffs_max, min_max = "max"):
    if len(coeffs_min) != len(coeffs_max):
        print("coeffs_min not same length as coeffs_max!")
        exit()
    if min_max == "max":
        return P_cot_delta(P, coeffs_max)
    elif min_max == "min":
        return P_cot_delta(P, coeffs_min)

    # results = np.zeros(2**len(coeffs_min))
    # for i in range(len(coeffs_min)):
    #     for j in range(2):
    #         coeffs_tmp = np.zeros(len(coeffs_min))
    #         for k in range(len(coeffs_min)):
    #             if j == 0:
    #                 coeffs_tmp[k] = coeffs_min[k]
    #             elif j == 1:
    #                 coeffs_tmp[k] = coeffs_max[k]
    #         results[2*i+j] = P_cot_delta(P, coeffs_tmp)
    # if min_max == "max":
    #     return max(results)
    # elif min_max == "min":
    #     return min(results)

def plot_P_cot_data(data, color, label):
    plt.errorbar(x=data[0], xerr = data[1], y=data[2], yerr=data[3], color = color, ls = "", capsize=5, markersize=10, label = label)

def plot_P_cot_fit(coeffs, coeffs_min, coeffs_max, color):
    size = 100
    Parr = np.linspace(0,1,size)
    P2arr = np.zeros(size)
    P_cot_arr = np.zeros(size)
    P_cot_arr_p = np.zeros(size)
    P_cot_arr_m = np.zeros(size)

    for i in range(len(P_cot_arr)):
        P2arr[i] = Parr[i]**2
        P_cot_arr[i] = P_cot_delta(Parr[i], coeffs)
        P_cot_arr_p[i] = P_cot_delta_min_max(Parr[i], coeffs_min, coeffs_max, min_max="max")
        P_cot_arr_m[i] = P_cot_delta_min_max(Parr[i], coeffs_min, coeffs_max, min_max="min")
    
    plt.plot(P2arr, P_cot_arr, color = color)
    plt.fill_between(x=P2arr, y1=P_cot_arr_m, y2=P_cot_arr_p, alpha = 0.3, color = color)

def plot_one_P_cot(filename, color, order):
    beta = float(filename[35:40])
    m_1 = float(filename[43:49])
    # print(beta)
    # print(m_1)
    # exit()
    input_data = get_data_from_file(beta, m_1, limit = 0.9, order=order)
    # input_data = get_data_from_file(filename, order)
    if input_data != None:
        data, coeffs, coeffs_p, coeffs_m, mass_Goldstone, beta, m = input_data
        plot_P_cot_data(data, color, "b69m09")
        plot_P_cot_fit(coeffs, coeffs_m, coeffs_p, color = color)
        plt.grid()
        plt.xlim((0,1.2))
        plt.ylim((-4,1))
        plt.xlabel("$P^{\prime 2}$")
        plt.ylabel("$P^{\prime}\cot(\delta)$")
        plt.title("beta = %1.3f, m1/2 = %1.3f, order = %i"%(beta, m, order))
        plt.savefig("plots/P_cot_PS_b%1.3f_m%1.3f_o%i.pdf"%(beta, m, order))
        plt.show()
        plt.clf()

def plot_all_P_cot():
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders     = [2,0,0,2,2,4,2,2,4,0,0,0]
    # orders_alt = [2,0,0,0,2,2,2,2,2,0,0,0]
    # "phase_shift_fit_P_cot_PS_SU(3)_beta5.400_m1-0.890_lim_0.9"
    for i, file in zip(range(len(filenames)), filenames):
        plot_one_P_cot(file, color = color_arr[i], order = orders[i])
        # for j in [0,2,4,6,8]:
        #     plot_one_P_cot(file, color = color_arr[i], order = j)
    # plt.grid()
    # plt.savefig("plots/P_cot_PS_all.pdf", bbox_inches = "tight")
    # plt.xscale("log")
    # plt.xlim([1e-3, 1])
    # plt.savefig("plots/P_cot_PS_all_log.pdf", bbox_inches = "tight")
    # plt.show()

def E_dispersion_relation(P, mass_Goldstone):       # care for correct use of lattice disp rel, takes P = P/mass_Goldstone, return E = E/mass_Goldstone
    return 2*np.sqrt(1+P**2)
    # return 2*np.arccosh(np.cosh(mass_Goldstone)+2*(np.sin(P*mass_Goldstone/2))**2)/mass_Goldstone

def sigma_of_P(P, P_cot, mass_Goldstone):                                              # P - momentum, coeffs contains coeffs of expansions
    # E = 2*np.arccosh(np.cosh(mass_Goldstone)+2*np.sin(P*mass_Goldstone/2)**2)/mass_Goldstone     # care for correct use of lattice disp rel
    E = E_dispersion_relation(P, mass_Goldstone)
    return 16*np.pi/(E**2*(1+(P_cot/P)**2))

def sigma(P, coeffs, mass_Goldstone):                                                   # takes real unitless P = P/mass_Goldstone
    return sigma_of_P(P, P_cot_delta(P, coeffs), mass_Goldstone)

def sigma_min_max(P, coeffs_min, coeffs_max, mass, min_max = "max"):
    return sigma_of_P(P, P_cot_delta_min_max(P, coeffs_min, coeffs_max, min_max=min_max), mass)

def P_of_s(s, mass_Goldstone):                              # care to use lattice disp correctly
    return 2*np.arcsin(np.sqrt((np.cosh(np.sqrt(s/4))-np.cosh(mass_Goldstone))/2))/mass_Goldstone

def s_of_P(P, mass_Goldstone):
    return E_dispersion_relation(P, mass_Goldstone)**2
    # return (2*np.arccosh(np.cosh(mass_Goldstone)+2*(np.sin(P*mass_Goldstone/2))**2)/mass_Goldstone)**2

def P_of_v(v, mass_Goldstone):                              # all real unitless
    return v*mass_Goldstone/np.sqrt(1-v*v)

def v_of_P(P, mass_Goldstone):
    return 1/np.sqrt(1+(mass_Goldstone/P)**2)

def plot_one_sigma(filename, color, order = 2, xaxis = "P"):
    data, coeffs, coeffs_p, coeffs_m, mass_Goldstone, beta, m = get_data_from_file(filename, order = order)
    num = 1000
    # P_max = P_of_s((4*mass_Goldstone)**2, mass_Goldstone)
    # Parr = np.linspace(0,P_max,num)
    Parr = np.linspace(0,3,num)                                                         # P here = P/mass_Goldstone
    sarr = np.linspace(0,10,num)
    varr = np.linspace(0,1,num)
    for i in range(len(Parr)):
        sarr[i] = s_of_P(Parr[i], mass_Goldstone)
        varr[i] = v_of_P(Parr[i], mass_Goldstone)
        # Parr[i] = P_of_v(varr[i], mass_Goldstone)
        # sarr[i] = s_of_P(P_of_v(varr[i], mass_Goldstone), mass_Goldstone)

    sigma_arr = np.zeros(num)
    sigma_arr_p = np.zeros(num)
    sigma_arr_m = np.zeros(num)

    for i in range(len(Parr)):
        sigma_arr[i] = sigma(Parr[i], coeffs, mass_Goldstone)
        sigma_arr_p[i] = sigma_min_max(Parr[i], coeffs_m, coeffs_p, mass_Goldstone, min_max="max")
        sigma_arr_m[i] = sigma_min_max(Parr[i], coeffs_m, coeffs_p, mass_Goldstone, min_max="min")

    if xaxis == "P":
        plt.plot(Parr, sigma_arr, color = color, label = "b%1.3f m%1.3f"%(beta, m))
        # plt.fill_between(Parr, sigma_arr_p, sigma_arr_m, color = color, alpha = 0.3)
    elif xaxis == "v":
        plt.plot(varr, sigma_arr, color = color, label = "b%1.3f m%1.3f"%(beta, m))
        # plt.fill_between(varr, sigma_arr_p, sigma_arr_m, color = color, alpha = 0.3)
    elif xaxis == "s":
        plt.plot(sarr, sigma_arr, color = color, label = "b%1.3f m%1.3f"%(beta, m))
        # plt.fill_between(sarr, sigma_arr_p, sigma_arr_m, color = color, alpha = 0.3)
    # plt.show()

def plot_all_sigma(xaxis = "P"):
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders     = [2,0,0,2,2,4,2,2,4,0,0,0]
    # orders_alt = [2,0,0,0,2,2,2,2,2,0,0,0]    
    # "phase_shift_fit_P_cot_PS_SU(3)_beta5.400_m1-0.890_lim_0.9"
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

def intfunc(v, vmean, coeffs, mass):
    return v*sigma(P_of_v(v, mass)/mass, coeffs, mass)*MB_dist(v, vmean)

    # tmp = v*sigma(P_of_v(v, mass)/mass, coeffs, mass)*MB_dist(v, vmean)
    # if tmp > 10:
    #     if v > 0.95:
    #         print()
    #         print(tmp)
    #         print(v)
    #         print(sigma(P_of_v(v, mass)/mass, coeffs, mass))
    #         print(MB_dist(v, vmean))
    #         print()
    # return tmp

def sigma_v_of_vmean(vmean, coeffs, mass):
    return scpint.quad(func=intfunc, a=0.0001, b=0.99999, args=(vmean,coeffs,mass))[0]    

def plot_one_sigma_v_vs_v(filename, color, order = 2, massscale = 0):
    data, coeffs, coeffs_p, coeffs_m, mass_Goldstone, beta, m = get_data_from_file(filename, order = order)
    # /(massscale**3*4.578*1e-6)*lightspeed
    num = 100
    vmeanarr = np.linspace(0, 1, num)
    sigma_v_arr = np.zeros(num)

    for i in range(len(vmeanarr)):
        sigma_v_arr[i] = sigma_v_of_vmean(vmeanarr[i], coeffs, mass_Goldstone)
        # tmp = sigma_v_of_vmean(vmeanarr[i], coeffs, mass_Goldstone)
        # if tmp < 1e3:
        #     sigma_v_arr[i] = tmp

    if massscale != 0:
        for i in range(len(vmeanarr)):
            vmeanarr[i] = vmeanarr[i] * lightspeed
            sigma_v_arr[i] = sigma_v_arr[i]/(massscale**3*4.578*1e-6)*lightspeed


    plt.plot(vmeanarr, sigma_v_arr, color = color, label = "b%1.3f m%1.3f"%(beta, m))
    # plt.show()

def plot_all_sigma_v(massscale = 0):
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders     = [2,0,0,2,2,4,2,2,4,0,0,0]
    # orders_alt = [2,0,0,0,2,2,2,2,2,0,0,0]    
    # "phase_shift_fit_P_cot_PS_SU(3)_beta5.400_m1-0.890_lim_0.9"
    for i, file in zip(range(len(filenames)), filenames):
        plot_one_sigma_v_vs_v(file, color = color_arr[i], order = orders[i], massscale=massscale)
    plt.grid()
    if massscale == 0:
        plt.xlabel("$<v>$")
        plt.ylabel("$<\sigma v>$/m")
    else: 
        HALOdata = np.transpose(np.genfromtxt("input/DMhalodata.dat"))
        plt.scatter(HALOdata[0], HALOdata[1], color = "grey", label = "Halo data")
        plt.xlabel("$<v>$ in km/s")
        plt.ylabel("$<\sigma v>$/m in cm^2/g km/s")
    # plt.legend()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig("plots/sigma_v_vmean_mDM_%i.pdf"%int(massscale), bbox_inches = "tight")
    # plt.show()
    plt.title("$m_{DM}$ = %i MeV"%int(massscale))
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("plots/sigma_v_vmean_log_mDM_%i.pdf"%int(massscale), bbox_inches = "tight")
    plt.show()

def plot_integrand(filename, color, order = 2):
    data, coeffs, coeffs_p, coeffs_m, mass_Goldstone, beta, m = get_data_from_file(filename, order=order)

    print(beta, m)

    num1 = 1000
    num2 = 5
    varr = np.linspace(0,1,num1)
    vmeanarr = np.linspace(0.1,0.9,num2)
    MB_dist_arr = np.zeros((num2,num1))
    intfunc_arr = np.zeros((num2,num1))

    v_min = 0.1
    v_max = 0.9

    for i in range(num2):
        for j in range(num1):
            # MB_dist_arr[i][j] = MB_dist(varr[j], vmeanarr[i])
            intfunc_arr[i][j] = intfunc(varr[j], vmeanarr[i], coeffs, mass_Goldstone)
        # plt.plot(varr, MB_dist_arr[i], label = "MB (v=%1.2f)"%vmeanarr[i], ls = "--")
        plt.plot(varr, intfunc_arr[i], label = "int_func (v=%1.2f)"%vmeanarr[i])

    plt.title("b%1.3f m%1.3f"%(beta, m))
    plt.grid()
    plt.legend()
    plt.savefig("plots/integrands_b_%1.3f m_%1.3f"%(beta, m)+".pdf", bbox_inches = "tight")
    # plt.show()
    plt.clf()

def plot_all_integrands():
    filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    orders     = [2,0,0,2,2,4,2,2,4,0,0,0]
    # orders_alt = [2,0,0,0,2,2,2,2,2,0,0,0]    
    # "phase_shift_fit_P_cot_PS_SU(3)_beta5.400_m1-0.890_lim_0.9"
    for i, file in zip(range(len(filenames)), filenames):
        plot_integrand(file, color = color_arr[i], order = orders[i])

    # plt.legend()
    # plt.savefig("plots/integrands.pdf")
    # plt.show()





def main():
    # for data in get_data_from_file("phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", 2):
    #     print(data)
    #     print()
    plot_all_P_cot()
    # plot_all_sigma(xaxis="P")
    plot_all_sigma(xaxis="v")
    # plot_all_sigma(xaxis="s")
    # plot_all_integrands()
    # plot_all_sigma_v(10)
    # plot_all_sigma_v(1000)
    plot_all_sigma_v(100)

    # filenames = ["phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.920_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.905_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.890_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.910_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.780_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.850_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.794_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.835_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.900_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta6.900_m1-0.870_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.200_m1-0.750_lim_0.9", "phase_shift_fit_P_cot_PS_SP(4)_beta7.050_m1-0.800_lim_0.9"]
    # orders     = [2,0,0,2,2,4,2,2,4,0,0,0]
    # orders_alt = [2,0,0,0,2,2,2,2,2,0,0,0]




if __name__ == "__main__":
    main()

















# def plot_sigma():
#     sigma = []
#     sigma_p = []
#     sigma_m = []
#     P = np.linspace(0,2,50)
#     # v = []
#     # E = []
#     # s = []
#     m_pi_rho_arr = []
#     beta_arr = []
#     m_arr = []
#     all_meas = errcl.measurement("all_SP4")
#     all_meas.read_from_HDF()
#     print("all meas read")

#     ind = 0

#     for result_name in all_meas.result_names:
#         search_str = "phase_shift_fit_P_cot_PS_SP(4)"
#         if result_name[:len(search_str)] == search_str:
#             search_str = "a2"
#             sigma.append([])
#             sigma_p.append([])
#             sigma_m.append([])
#             if result_name[len(result_name)-len(search_str):] == search_str:
#                 search_str = "lim_0/"                                               # which limit for P do you want?
#                 if search_str in result_name:
#                     info_str = result_name[:len(result_name)-2]
#                     beta = all_meas.infos[info_str+"beta"]
#                     beta_arr.append(beta)
#                     m = all_meas.infos[info_str+"m_1"]
#                     m_arr.append(m)
#                     m_pi_rho = all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].median[0]
#                     m_pi_rho_arr.append(m_pi_rho)

#                     mass = all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/mass_Goldstone"%(beta,m,m)].median[0]
#                     coeffs = []
#                     coeffs.append(all_meas.results[result_name[:len(result_name)-2]+"a2"].median[0])
#                     coeffs.append(all_meas.results[result_name[:len(result_name)-2]+"b2"].median[0])
#                     for i in range(len(P)):
#                         sigma[ind].append(sigma_of_P(P[i], coeffs, mass))
#                     plt.plot(P, sigma[ind])
#                     print(len(sigma[ind]))

#                     ind += 1
#                     print(result_name)
#     plt.show()




    # data = {}
    # beta_m_arr = []
    # all_meas = errcl.measurement("all_SP4")
    # all_meas.read_from_HDF()
    # print("all meas read")

    # ind = 0

    # for result_name in all_meas.result_names:
    #     search_str = "phase_shift_fit_P_cot_PS_SP(4)"
    #     if result_name[len(result_name)-len(search_str):] == search_str:
    #         search_str = "lim_0/"                                               # which limit for P do you want?
    #         if search_str in result_name:
    #             info_str = result_name[:len(result_name)-2]
    #             beta = all_meas.infos[info_str+"beta"]
    #             m = all_meas.infos[info_str+"m_1"]
    #             beta_m_str = "beta%1.3f_m%1.3f"%(beta,m)
    #             beta_m_arr.append(beta_m_str)
    #             data[beta_m_str] = {}
    #             data[beta_m_str]["m"] = m
    #             data[beta_m_str]["beta"] = beta
    #             for i in range(10):
    #                 if result_name[len(result_name)-i] == "/":
    #                     data[beta_m_str][result_name[len(result_name)-i:]] = all_meas.results[result_name].median[0]
    #                     data[beta_m_str][result_name[len(result_name)-i:]+"ep"] = all_meas.results[result_name].ep[0]
    #                     data[beta_m_str][result_name[len(result_name)-i:]+"em"] = all_meas.results[result_name].em[0]





    #             m_pi_rho = all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].median[0]
    #             m_pi_rho_arr.append(m_pi_rho)
    #             m_pi = all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/mass_Goldstone"%(beta,m,m)].median[0]
    #             m_pi_arr.append(m_pi)
    #             coeffs = []
    #             coeffs.append(all_meas.results[result_name[:len(result_name)-2]+"a2"].median[0])
    #             coeffs.append(all_meas.results[result_name[:len(result_name)-2]+"b2"].median[0])