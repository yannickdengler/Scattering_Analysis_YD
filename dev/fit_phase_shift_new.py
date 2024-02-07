import numpy as np
import error_classes as errcl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
# import os
import h5py
# from scipy.integrate import quad as integrate
import random
import math


def fit_func_const(P2, a):
    return a + 0*P2
def fit_func_p2(P2, a, b):
    return a + b*P2
def fit_func_p4(P2, a, b, c):
    return a + b*P2 + c*P2*P2
def fit_func_p6(P2, a, b, c, d):
    return a + b*P2 + c*P2*P2 + d*P2*P2*P2

def fit_func(order):
    if order == 0:
        return fit_func_const
    elif order == 2:
        return fit_func_p2
    elif order == 4:
        return fit_func_p4
    elif order == 6:
        return fit_func_p6
    else:
        print("Order given to fit_func is not correct!")
        exit()

def get_fit_params(P_2, P_cot_PS, num_ensembles):
    result = {}
    if num_ensembles > 1:
        popt, pcov = curve_fit(fit_func(0), P_2, P_cot_PS)
        a0 = popt[0]
        result["a0"] = [a0,]
        result["a_0_0"] = [1./a0,]
        if num_ensembles > 2:
            popt, pcov = curve_fit(fit_func(2), P_2, P_cot_PS)
            a2, b2 = popt
            result["a2"], result["b2"] = [[a2,],[b2,]]
            result["a_0_2"] = [1./a2,]
            # xarr = np.linspace(0.5*min(P_2), 2*max(P_2))
            # yarr = np.zeros(len(xarr))
            # for i in range(len(xarr)):
            #     yarr[i] = fit_func(2)(xarr[i], a2,b2)
            # plt.xlim([0,0.9])
            # plt.ylim([-5,5])
            # plt.scatter(P_2, P_cot_PS)
            # plt.plot(xarr, yarr)
            # plt.show()
            if num_ensembles > 3:
                popt, pcov = curve_fit(fit_func(4), P_2, P_cot_PS)
                a4, b4, c4 = popt
                result["a4"], result["b4"], result["c4"] = [[a4,],[b4,],[c4,]]
                result["a_0_4"] = [1./a4,]
                # xarr = np.linspace(0.5*min(P_2), 2*max(P_2))
                # yarr = np.zeros(len(xarr))
                # for i in range(len(xarr)):
                #     yarr[i] = fit_func(4)(xarr[i], a4,b4,c4)
                # plt.xlim([0,0.9])
                # plt.ylim([-5,5])
                # plt.scatter(P_2, P_cot_PS)
                # plt.plot(xarr, yarr)
                # plt.show()
                if num_ensembles > 4:
                    popt, pcov = curve_fit(fit_func(6), P_2, P_cot_PS)
                    a6, b6, c6, d6 = popt
                    result["a6"], result["b6"], result["c6"], result["d6"] = [[a6,],[b6,],[c6,],[d6,]]
                    result["a_0_6"] = [1./a6,]
    return result

def get_data_from_file(filelist, limit = 0, levels = [0,]):
    P_2 = []
    P_cot_PS = []
    filenames = []
    montecarlo_times = []
    N_Ls = []
    mass_Goldstone = np.zeros(3)
    beta = None
    m_1 = None
    P_2_median = []
    P_2_ep = []
    P_2_em = []
    P_cot_PS_median = []
    P_cot_PS_ep = []
    P_cot_PS_em = []
    for resultfile in filelist:
        if int(resultfile[25]) in levels:
            phase_shift = errcl.measurement(resultfile)
            if phase_shift.file_exists():
                phase_shift.read_from_HDF()
                if phase_shift.results["P_2_prime"].median[0] < limit or limit == 0:
                    N_Ls.append(phase_shift.results["N_L"].median[0])
                    mass_Goldstone[0] = phase_shift.results["mass_Goldstone"].median[0]
                    mass_Goldstone[1] = phase_shift.results["mass_Goldstone"].ep[0]
                    mass_Goldstone[2] = phase_shift.results["mass_Goldstone"].em[0]
                    infos = phase_shift.infos
                    # print(resultfile)
                    # print(infos["N_mont"])
                    infos["levels_fit"] = levels
                    infos["limit_fit"] = limit
                    beta = infos["beta"]
                    m_1 = infos["m_1"]
                    P_cot_PS_tmp = np.swapaxes(phase_shift.results["P_cot_PS_prime"].sample,0,1)[0]                             # check if you want prime or not
                    P_2_tmp = np.swapaxes(phase_shift.results["P_2_prime"].sample,0,1)[0]
                    if not(any(np.isnan(P_cot_PS_tmp)) or any(np.isinf(P_cot_PS_tmp)) or any(np.isnan(P_2_tmp)) or any(np.isinf(P_2_tmp))):
                        filenames.append(str(resultfile))
                        montecarlo_times.append(infos["N_mont"])
                        P_cot_PS.append(np.swapaxes(phase_shift.results["P_cot_PS_prime"].sample,0,1)[0])
                        P_2.append(np.swapaxes(phase_shift.results["P_2_prime"].sample,0,1)[0])
                        P_2_median.append(phase_shift.results["P_2_prime"].median[0])
                        P_2_ep.append(phase_shift.results["P_2_prime"].ep[0])
                        P_2_em.append(phase_shift.results["P_2_prime"].em[0])
                        P_cot_PS_median.append(phase_shift.results["P_cot_PS_prime"].median[0])
                        P_cot_PS_ep.append(phase_shift.results["P_cot_PS_prime"].ep[0])
                        P_cot_PS_em.append(phase_shift.results["P_cot_PS_prime"].em[0])
                    else:
                        print("Nans of infs in: "+resultfile)
    plot_data = [P_2_median, P_2_ep, P_2_em, P_cot_PS_median, P_cot_PS_ep, P_cot_PS_em]
    # print(montecarlo_times)
    return P_2, P_cot_PS, filenames, montecarlo_times, beta, m_1, mass_Goldstone, N_Ls, plot_data

def weight_by_montecarlotimes(P_2, P_cot_PS, montecarlotimes, num_reweight = 2000):
    num_ensembles = len(P_2)
    P_2_weigthed = np.zeros(num_reweight)
    P_cot_PS_weigthed = np.zeros(num_reweight)
    sum_montecarlotimes = sum(montecarlotimes)
    montecarlotime_arr = np.zeros(len(montecarlotimes))
    for i in range(len(montecarlotime_arr)):
        if i == 0:
            montecarlotime_arr[i] = montecarlotimes[i]
        else:
            montecarlotime_arr[i] = montecarlotime_arr[i-1] + montecarlotimes[i]

    for i in range(num_reweight):
        randint1 = random.randint(0, sum_montecarlotimes-1)
        for j in range(len(montecarlotime_arr)):
            if randint1 < montecarlotime_arr[j]:
                # print(j)
                randint2 = random.randint(0, len(P_2[j])-1)
                P_2_weigthed[i] = P_2[j][randint2]
                P_cot_PS_weigthed[i] = P_cot_PS[j][randint2]
    return P_2_weigthed, P_cot_PS_weigthed

def get_plus_minus(P_2, P_cot_PS, num_pm = 100):
    num_tot = len(P_2)
    if len(P_cot_PS) != num_tot:
        print("len of P_2 and P_cot_PS not same in get_plus_minus")
        exit()
    P_2_sorted = [x for _, x in sorted(zip(P_cot_PS,P_2))]
    P_cot_PS_sorted = sorted(P_cot_PS)
    percentage_std = 0.682689
    ind_mean_high = math.ceil(num_tot/2+num_pm/2)
    ind_mean_low = math.floor(num_tot/2-num_pm/2)
    ind_p_high = math.ceil(num_tot*percentage_std+num_pm/2)
    ind_p_low = math.floor(num_tot*percentage_std-num_pm/2)
    ind_m_high = math.ceil(num_tot*(1-percentage_std)+num_pm/2)
    ind_m_low = math.floor(num_tot*(1-percentage_std)-num_pm/2)
    P_2_mean = P_2_sorted[ind_mean_low:ind_mean_high]
    P_2_m = P_2_sorted[ind_m_low:ind_m_high]
    P_2_p = P_2_sorted[ind_p_low:ind_p_high]
    P_cot_PS_mean = P_cot_PS_sorted[ind_mean_low:ind_mean_high]
    P_cot_PS_m = P_cot_PS_sorted[ind_m_low:ind_m_high]
    P_cot_PS_p = P_cot_PS_sorted[ind_p_low:ind_p_high]
    return P_2_mean,P_2_m, P_2_p,P_cot_PS_mean, P_cot_PS_m, P_cot_PS_p

def get_median_and_err_from_P_2_P_cot_PS(P_2, P_cot_PS):
    percentage_std = 0.682689
    P_2_mean = []
    P_2_ep = []
    P_2_em = []
    P_cot_PS_mean = []
    P_cot_PS_ep = []
    P_cot_PS_em = []
    for i in range(len(P_2)):
        low_ind = math.ceil(num*(1-percentage_std)/2)
        high_ind = math.floor(num*(1+percentage_std)/2)
        P_2_mean.append(P_2[i][len(P_2[i]//2)])
        P_2_ep.append(P_2[i][len(P_2[i]//2)])
        P_2_em.append(P_2[i][len(P_2[i]//2)])


def fit_pm_P_cot_PS(filelist, limit = 0, levels = [0,]):
    result = {}
    P_2, P_cot_PS, filenames, montecarlo_times, beta, m_1, mass_Goldstone, N_Ls, plot_data = get_data_from_file(filelist=filelist, limit = limit, levels = levels)

    result["N_Ls"] = N_Ls
    result["mass_Goldstone_median"] = [mass_Goldstone[0],]
    result["mass_Goldstone_ep"] = [mass_Goldstone[1],]
    result["mass_Goldstone_em"] = [mass_Goldstone[2],]
    result["plot_data"] = plot_data
    result["beta"] = [beta,]
    result["m_1"] = [m_1,]
    print(beta, m_1)
    num_ensembles = len(P_2)

    P_2_median = []
    P_2_m = []
    P_2_p = []
    P_cot_PS_median = []
    P_cot_PS_m = []
    P_cot_PS_p = []
    for i in range(num_ensembles):
        P_2_median_tmp, P_2_m_tmp, P_2_p_tmp, P_cot_PS_median_tmp, P_cot_PS_m_tmp, P_cot_PS_p_tmp = get_plus_minus(P_2[i], P_cot_PS[i], num_pm=50)
        P_2_median.append(P_2_median_tmp)
        P_2_m.append(P_2_m_tmp)
        P_2_p.append(P_2_p_tmp)
        P_cot_PS_median.append(P_cot_PS_median_tmp)
        P_cot_PS_m.append(P_cot_PS_m_tmp)
        P_cot_PS_p.append(P_cot_PS_p_tmp)
    
    P_2_orig = []
    P_cot_PS_orig = []
    P_2_rw = []
    P_cot_PS_rw = []
    P_2_m_rw = []
    P_2_p_rw = []
    P_cot_PS_m_rw = []
    P_cot_PS_p_rw = []

    for i in range(num_ensembles):                                    # THIS WAS FOR TESTING IF REWEIGHTING WORKS!!
        # for j in range(len(P_2[i])):
            # P_2_rw.append(P_2[i][j])
            # P_cot_PS_rw.append(P_cot_PS[i][j])
        for j in range(len(P_2[i])):
            P_2_orig.append(P_2[i][j])
            P_cot_PS_orig.append(P_cot_PS[i][j])
        for j in range(len(P_2_median[i])):
            P_2_rw.append(P_2_median[i][j])
            P_cot_PS_rw.append(P_cot_PS_median[i][j])
        for j in range(len(P_2_m[i])):
            P_2_m_rw.append(P_2_m[i][j])
            P_cot_PS_m_rw.append(P_cot_PS_m[i][j])
        for j in range(len(P_2_p[i])):
            P_2_p_rw.append(P_2_p[i][j])
            P_cot_PS_p_rw.append(P_cot_PS_p[i][j])

    # for i in range(num_ensembles):
    #     P_2_rw, P_cot_PS_rw = weight_by_montecarlotimes(P_2, P_cot_PS, montecarlo_times, num_reweight=2000)
    #     P_2_m_rw, P_cot_PS_m_rw = weight_by_montecarlotimes(P_2_m, P_cot_PS_m, montecarlo_times, num_reweight=1000)
    #     P_2_p_rw, P_cot_PS_p_rw = weight_by_montecarlotimes(P_2_p, P_cot_PS_p, montecarlo_times, num_reweight=1000)

    # print(montecarlo_times)

    results_all = get_fit_params(P_2_orig, P_cot_PS_orig, num_ensembles)
    results_mean = get_fit_params(P_2_rw, P_cot_PS_rw, num_ensembles)
    results_m = get_fit_params(P_2_m_rw, P_cot_PS_m_rw, num_ensembles)
    results_p = get_fit_params(P_2_p_rw, P_cot_PS_p_rw, num_ensembles)
    plt.show()
    plt.clf()
    print("show")

    for key, val in results_mean.items():
        result[key+"_mean"] = val
    for key, val in results_m.items():
        result[key+"_m"] = val
    for key, val in results_p.items():
        result[key+"_p"] = val
    with h5py.File("/home/dengler_yannick/Documents/Scattering_Analysis_YD/output/phase_shift_fit_results/phase_shift_fit_results_b%1.3f_m%1.3f_lim%1.3f"%(beta,m_1, limit),"w") as f:
        for key, val in result.items():
            f.create_dataset(key, data = val)
    # with h5py.File("/home/dengler_yannick/Documents/Scattering_Analysis_YD/output/phase_shift_fit_results/phase_shift_fit_results_b%1.3f_m%1.3f_lim%1.3f"%(beta,m_1, limit),"r") as f:
    #     for key in f.keys():
    #         print(key, f[key][:])

def create_all_filenames():
    beta_m_str_list = []
    PATH = "output/result_files/"
    temp = "phase_shift_Fabian"
    filelist = os.listdir(PATH)
    resultfile_list = []
    num = len(temp)
    for file in filelist:
        length = len(file)
        if file[:num] == temp:
            print(file)
            beta_m_str = file[47:74]
            if not(beta_m_str in beta_m_str_list):
                beta_m_str_list.append(beta_m_str)
                resultfile_list.append([])
                resultfile_list[beta_m_str_list.index(beta_m_str)].append(file[:length-5])
            else:
                resultfile_list[beta_m_str_list.index(beta_m_str)].append(file[:length-5])
    for listi in resultfile_list:
        print(listi)
    with open("input/filenames_phase_shift_fit_all", "w") as filestream:
        for filelist in resultfile_list:
            for file in filelist:
                filestream.write(file+"\t")
            filestream.write("\n")

def main():
    with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_fit", "r") as f:
        RL = f.readlines()
        for filelist in RL:
            # print(filelist.split())
            # print()
            # print()
            fit_pm_P_cot_PS(filelist.split(), limit = 0.9)

if __name__ == "__main__":
    # create_all_filenames()
    main()

