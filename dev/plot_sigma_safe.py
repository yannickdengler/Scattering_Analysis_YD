import os
import matplotlib.pyplot as plt
import read_HDF5_logfile as HDF_log
import error_classes as errcl 
import numpy as np
import matplotlib
from matplotlib import cm
import math

from scipy.interpolate import interp1d as inter

color_arr = ["blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick"]

font = {'size'   : 16}
matplotlib.rc('font', **font)

def fit_P_cot_PS(filelist, x, y, error, color = "red"):
    bm_str = filelist[0][41:65]
    meas = errcl.measurement("phase_shift_fit_P_cot_PS_"+bm_str)
    if meas.file_exists():
        meas.read_from_HDF()
        fit_result = meas.results
        num = len(meas.results[x].median)
        # xarr = np.linspace(min(meas.results[x].median)*0.9,max(meas.results[x].median)*1.1, 500)
        xarr = np.linspace(0,2, 500)
        yarr_2 = fit_func_p2(xarr, fit_result["a2"].median[0], fit_result["b2"].median[0])
        plt.plot(xarr, yarr_2, color = color, label = "O(P2)")
        if error:
            yarr_2m = fit_func_p2(xarr, fit_result["a2"].median_m[0], fit_result["b2"].median[0])  
            yarr_2p = fit_func_p2(xarr, fit_result["a2"].median_p[0], fit_result["b2"].median[0])  
            plt.fill_between(x=xarr, y1=yarr_2m, y2=yarr_2p, color = color, alpha = 0.1)
        # if num > 2:
        #     yarr_4 = fit_func_p4(xarr, fit_result["a4"].median[0], fit_result["b4"].median[0], fit_result["c4"].median[0])
        #     plt.plot(xarr, yarr_4, color = color, ls = "--", label = "O(P4)")
        #     if error:
        #         yarr_4m = fit_func_p4(xarr, fit_result["a4"].median_m[0], fit_result["b4"].median[0], fit_result["c4"].median[0])  
        #         yarr_4p = fit_func_p4(xarr, fit_result["a4"].median_p[0], fit_result["b4"].median[0], fit_result["c4"].median[0])  
        #         plt.fill_between(x=xarr, y1=yarr_4m, y2=yarr_4p, color = color, alpha = 0.1)
        print(bm_str)
        print("a_0 (order 1) = %1.2e + %1.2e - %1.2e"%(fit_result["a_0_2"].median[0],fit_result["a_0_2"].ep[0],fit_result["a_0_2"].em[0]))
        if fit_result["a_0_2"].median[0] + fit_result["a_0_2"].ep[0] > 0:
            print("consistent with 0 in 1 sigma")
        if num > 2:
            print("a_0 (order 2) = %1.2e + %1.2e - %1.2e"%(fit_result["a_0_4"].median[0],fit_result["a_0_4"].ep[0],fit_result["a_0_4"].em[0]))
            if fit_result["a_0_4"].median[0] + fit_result["a_0_4"].ep[0] > 0:
                print("consistent with 0 in 1 sigma")
        print("")
        
def marker_beta(beta):
    if beta == 6.9:
        return "^"
    if beta == 7.05:
        return "s"
    if beta == 7.2:
        return "o"
        
def color_beta(beta):
    if beta == 6.9:
        return "#4363d8"
    if beta == 7.05:
        return "#800000"
    if beta == 7.2:
        return "#f032e6"

def plot_a_0_vs_m_f_pi(squared):
    plt.figure(figsize=(8,4.8))
    def chipt(x):
        return -x/32
    data = np.genfromtxt("output/plot_files/chipt.dat")
    # data = np.genfromtxt("output/plot_files/chipt_only_stat.dat")
    for i in range(len(data)):
        if squared:
            plt.errorbar(x=data[i][5],xerr=data[i][6], y=[data[i][0],],yerr=[[data[i][1],],[data[i][2],]], marker = marker_beta(data[i][7]), ls = "", capsize=5, markersize=10, color = color_beta(data[i][7]))
        else:
            plt.errorbar(x=data[i][3],xerr=data[i][4], y=[data[i][0],],yerr=[[data[i][1],],[data[i][2],]], marker = marker_beta(data[i][7]), ls = "", capsize=5, markersize=10, color = color_beta(data[i][7]))
    plt.fill_between(x=(-100,100), y1=(-0.43381100719774743,-0.43381100719774743), y2=(-0.9218134335278696,-0.9218134335278696), color = "grey", alpha = 0.25)
    plt.axhline(-0.6444469758057019, color = "grey", alpha = 0.5)

    if squared:
        xarr = np.linspace(0,100)
        yarr = [chipt(x) for x in xarr]
    else:
        xarr = np.linspace(0,10)
        yarr = [chipt(x**2) for x in xarr]
    
    plt.plot(xarr, yarr, ls = "dashed", color = "green", label = "LO EFT")

    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(6.9), color = color_beta(6.9), label ="$\\beta=6.90$", s = 60)
    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(7.05), color = color_beta(7.05), label ="$\\beta=7.05$", s = 60)
    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(7.2), color = color_beta(7.2), label ="$\\beta=7.20$", s = 60)

    if squared:
        # plt.xlim([18,50])
        plt.xlim([0,50])
    else:
        # plt.xlim([4,7.5])
        plt.xlim([0,7.5])
    # plt.ylim([-1.02,-0.18])
    plt.ylim([-1.4,0])
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend()
    if squared:
        plt.xlabel("$\\frac{m_\pi^2}{f_\\pi^2}$")
    else:
        plt.xlabel("$\\frac{m_\pi}{f_\\pi}$")
    plt.ylabel("$a_0m_\pi$")
    plt.grid()
    if squared:
        # plt.savefig("plots/a_0_vs_m_f_pi_2.pdf", bbox_inches="tight")
        plt.savefig("plots/a_0_vs_m_f_pi_2_nozoom.pdf", bbox_inches="tight")
    else:
        # plt.savefig("plots/a_0_vs_m_f_pi.pdf", bbox_inches="tight")
        plt.savefig("plots/a_0_vs_m_f_pi_nozoom.pdf", bbox_inches="tight")
    plt.show()
    plt.clf()

def plot_a_0_vs_m_pi_m_rho():
    a_0_2 = []
    a_0_2_p = []
    a_0_2_m = []
    m_pi_rho_2 = []
    m_pi_rho_2_p = []
    m_pi_rho_2_m = []
    a_0_4 = []
    a_0_4_p = []
    a_0_4_m = []
    m_pi_rho_4 = []
    m_pi_rho_4_p = []
    m_pi_rho_4_m = []
    beta_arr = []
    m_arr = []
    all_meas = errcl.measurement("all_SP4")
    all_meas.read_from_HDF()

    a_0_total = []

    for result_name in all_meas.result_names:
        search_str = "phase_shift_fit_P_cot_PS_SP(4)"
        if result_name[:len(search_str)] == search_str:
            search_str = "a_0_2"
            if result_name[len(result_name)-len(search_str):] == search_str:
                print(result_name)
                a_0_2.append(all_meas.results[result_name].median)
                a_0_2_p.append(all_meas.results[result_name].ep)
                a_0_2_m.append(all_meas.results[result_name].em)
                beta = float(result_name[35:40])
                beta_arr.append(beta)
                m = float(result_name[43:49])
                m_arr.append(m)
                m_pi_rho = all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].median[0]
                for i in range(len(all_meas.results[result_name].sample)):
                    a_0_total.append(all_meas.results[result_name].sample[i][0])
                m_pi_rho_2.append(all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].median)
                m_pi_rho_2_p.append(all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].ep)
                m_pi_rho_2_m.append(all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].em)
                print("beta = %1.3f, m = %1.3f, m_pi_rho = %1.3f"%(beta, m, m_pi_rho))
            search_str = "a_0_4"
            if result_name[len(result_name)-len(search_str):] == search_str:
                a_0_4.append(all_meas.results[result_name].median)
                a_0_4_p.append(all_meas.results[result_name].ep)
                a_0_4_m.append(all_meas.results[result_name].em)
                beta = float(result_name[35:40])
                m = float(result_name[43:49])
                m_pi_rho_4.append(all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].median)
                m_pi_rho_4_p.append(all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].ep)
                m_pi_rho_4_m.append(all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].em)
    percentage_std = 0.682689
    num = len(a_0_total)                                                  # num_results
    tmp = np.sort(a_0_total, kind="mergesort")
    a_0_median = tmp[num//2]
    low_ind = math.ceil(num*(1-percentage_std)/2)
    high_ind = math.floor(num*(1+percentage_std)/2)
    a_0_ep = abs(tmp[high_ind]-a_0_median)
    a_0_em = abs(tmp[low_ind]-a_0_median)
    print("a_0_mean: ", a_0_median, " + ", a_0_ep, " - ", a_0_em)
    plt.axhline(a_0_median, color = "grey", alpha = 0.5)
    plt.fill_between(x=(0,3), y1=(a_0_median+a_0_ep,a_0_median+a_0_ep), y2=(a_0_median-a_0_em,a_0_median-a_0_em), color = "grey", alpha = 0.25)
    beta_used = []
    for i in range(len(a_0_2)):                
        plt.errorbar(x=m_pi_rho_2[i],xerr=[m_pi_rho_2_p[i],m_pi_rho_2_m[i]], y=a_0_2[i],yerr=[a_0_2_p[i],a_0_2_m[i]], marker = marker_beta(beta_arr[i]), ls = "", capsize=5, markersize=10, color = color_beta(beta_arr[i]))

        # if beta_arr[i] in beta_used:
        #     plt.errorbar(x=m_pi_rho_2[i],xerr=[m_pi_rho_2_p[i],m_pi_rho_2_m[i]], y=a_0_2[i],yerr=[a_0_2_p[i],a_0_2_m[i]], marker = marker_beta(beta_arr[i]), ls = "", capsize=5, markersize=10, color = color_beta(beta_arr[i]))
        # else:
        #     plt.errorbar(x=m_pi_rho_2[i],xerr=[m_pi_rho_2_p[i],m_pi_rho_2_m[i]], y=a_0_2[i],yerr=[a_0_2_p[i],a_0_2_m[i]], marker = marker_beta(beta_arr[i]), ls = "", capsize=5, markersize=10, color = color_beta(beta_arr[i]), label = "$\\beta$ = %1.3f"%beta_arr[i])
        #     beta_used.append(beta_arr[i])

    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(6.9), color = color_beta(6.9), label ="$\\beta=6.90$", s = 60)
    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(7.05), color = color_beta(7.05), label ="$\\beta=7.05$", s = 60)
    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(7.2), color = color_beta(7.2), label ="$\\beta=7.20$", s = 60)

    plt.xlim([0.63,0.93])
    plt.ylim([-1.02,-0.18])
    # handles, labels = plt.gca().get_legend_handles_labels()
    # order = [1,0,2]
    # plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
    plt.legend()
    plt.xlabel("$\\frac{m_\pi}{m_\\rho}$")
    plt.ylabel("$a_0m_\pi$")
    plt.grid()
    plt.savefig("plots/a_0_vs_m_pi.pdf", bbox_inches="tight")
    plt.show()
    plt.clf()

    # for result_name in all_meas.result_names:
    #     search_str = "phase_shift_fit_P_cot_PS_SP(4)"
    #     if result_name[:len(search_str)] == search_str:
    #         search_str = "a_0_2"
    #         if result_name[len(result_name)-len(search_str):] == search_str:
    #             beta = float(result_name[35:40])
    #             beta_arr.append(beta)
    #             m = float(result_name[43:49])
    #             m_arr.append(m)
    #             m_pi_rho = all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].median[0]
    #             m_pi_rho_2.append(m_pi_rho)
    #             # print(m_pi_rho)
    #             m_pi_rho_2_p.append(all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].ep[0])
    #             m_pi_rho_2_m.append(all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].em[0])
    #             a_0_2.append(all_meas.results[result_name].median[0]/m_pi_rho)
    #             a_0_2_p.append(all_meas.results[result_name].ep[0]/m_pi_rho)
    #             a_0_2_m.append(all_meas.results[result_name].em[0]/m_pi_rho)
    #             print("beta = %1.3f, m = %1.3f, m_pi_rho = %1.3f"%(beta, m, m_pi_rho))
    #         search_str = "a_0_4"
    #         if result_name[len(result_name)-len(search_str):] == search_str:
    #             beta = float(result_name[35:40])
    #             m = float(result_name[43:49])
    #             m_pi_rho_4.append(all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].median[0])
    #             m_pi_rho_4_p.append(all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].ep[0])
    #             m_pi_rho_4_m.append(all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].em[0])
    #             a_0_4.append(all_meas.results[result_name].median[0]/m_pi_rho)
    #             a_0_4_p.append(all_meas.results[result_name].ep[0]/m_pi_rho)
    #             a_0_4_m.append(all_meas.results[result_name].em[0]/m_pi_rho)

    # plt.errorbar(x=m_pi_rho_2,xerr=[m_pi_rho_2_p,m_pi_rho_2_m], y=a_0_2,yerr=[a_0_2_p,a_0_2_m], ls = "", capsize=5, markersize=10, color = "red", label = "O(2)", marker = "^")
    # plt.errorbar(x=m_pi_rho_4,xerr=[m_pi_rho_4_p,m_pi_rho_4_m], y=a_0_4,yerr=[a_0_4_p,a_0_4_m], ls = "", capsize=5, markersize=10, color = "blue", label = "O(4)", marker = "v") 
    # plt.grid()
    # plt.legend()
    # for i in range(len(a_0_2)):
    #     plt.text(x=m_pi_rho_2[i]+0.01, y = a_0_2[i]-0.1, s = "b=%1.3f\nm=%1.3f"%(beta_arr[i],m_arr[i]), fontsize=6)
    # plt.title("$a_0$ vs $\\frac{m_\pi}{m_\\rho}$ with O(2) and O(4)-fits")
    # plt.xlabel("$\\frac{m_\pi}{m_\\rho}$")
    # plt.ylabel("$a_0m_\\rho$")
    # plt.savefig("plots/a_0_vs_m_rho.pdf")
    # plt.show()

    # plt.errorbar(x=beta_arr, y=a_0_2,yerr=[a_0_2_p,a_0_2_m], ls = "", capsize=5, markersize=10, color = "green", label = "$a_0$", marker = "^")
    # plt.grid()
    # plt.legend()
    # for i in range(len(a_0_2)):
    #     plt.text(x=beta_arr[i]+0.02, y = a_0_2[i], s = "m=%1.3f"%(m_arr[i]), fontsize=6)
    # plt.title("$a_0$ vs $\\beta$ with O(2)-fits")
    # plt.xlabel("$\\beta$")
    # plt.ylabel("$a_0m_\pi$")
    # plt.savefig("plots/a_0_vs_beta.pdf")
    # plt.show()

    # plt.errorbar(x=m_arr, y=a_0_2,yerr=[a_0_2_p,a_0_2_m], ls = "", capsize=5, markersize=10, color = "purple", label =  "$a_0$", marker = "^")
    # plt.grid()
    # plt.legend()
    # for i in range(len(a_0_2)):
    #     plt.text(x=m_arr[i]+0.002, y = a_0_2[i]+0.03, s = "m=%1.3f"%(beta_arr[i]), fontsize=6)
    # plt.title("$a_0$ vs m (bare fermion mass) with O(2)-fits")
    # plt.xlabel("$m_{1/2}$ (bare fermion mass)")
    # plt.ylabel("$a_0m_\pi$")
    # plt.savefig("plots/a_0_vs_m.pdf")
    # plt.show()

def plot_P_cot_PS():
    P_cot_PS = []
    P_cot_PS_p = []
    P_cot_PS_m = []
    a2 = []
    b2 = []
    P_2 = []
    P_2_p = []
    P_2_m = []
    m_pi_rho_arr = []
    beta_arr = []
    m_arr = []
    all_meas = errcl.measurement("all_SP4")
    all_meas.read_from_HDF()

    for result_name in all_meas.result_names:
        search_str = "phase_shift_Fabian_level_0_Scattering_I2_SP(4)"
        # print(result_name)
        if result_name[:len(search_str)] == search_str:
            search_str = "P_cot_PS_prime"
            if result_name[len(result_name)-len(search_str):] == search_str:
                info_str = result_name[:len(result_name)-len(search_str)]
                P_cot_PS.append(all_meas.results[result_name].median)
                P_cot_PS_p.append(all_meas.results[result_name].ep)
                P_cot_PS_m.append(all_meas.results[result_name].em)
                P_2.append(all_meas.results[result_name[:len(result_name)-len(search_str)]+"P_2_prime"].median)
                P_2_p.append(all_meas.results[result_name[:len(result_name)-len(search_str)]+"P_2_prime"].ep)
                P_2_m.append(all_meas.results[result_name[:len(result_name)-len(search_str)]+"P_2_prime"].em)
                beta = all_meas.infos[info_str+"beta"]
                beta_arr.append(beta)
                m = all_meas.infos[info_str+"m_1"]
                m_arr.append(m)
                m_pi_rho = all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].median[0]
                m_pi_rho_arr.append(m_pi_rho)
                # print("beta = %1.3f, m = %1.3f, m_pi_rho = %1.3f"%(beta, m, m_pi_rho))
        search_str = "phase_shift_fit_P_cot_PS_SP(4)"
        if result_name[:len(search_str)] == search_str:
            search_str = "a2"
            if result_name[len(result_name)-len(search_str):] == search_str:
                info_str = result_name[:len(result_name)-len(search_str)]
                a2.append(all_meas.results[result_name].median)
                b2.append(all_meas.results[result_name[:len(result_name)-len(search_str)]+"b2"].median)


    norm = matplotlib.colors.Normalize(vmin=min(m_pi_rho_arr), vmax=max(m_pi_rho_arr), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    time_color = np.array([(mapper.to_rgba(v)) for v in m_pi_rho_arr])
    plt.colorbar(mappable=mapper, label = "$\\frac{m_\pi}{m_\\rho}$")
    beta_used = []
    for i in range(len(P_2)):   
        plt.errorbar(x=P_2[i], xerr=[P_2_m[i],P_2_p[i]], y=P_cot_PS[i], yerr=[P_cot_PS_m[i],P_cot_PS_p[i]], marker = marker_beta(beta_arr[i]), color = time_color[i], capsize=5, markersize=10)


    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(6.9), color = "grey", label ="$\\beta=6.90$", s = 60)
    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(7.05), color = "grey", label ="$\\beta=7.05$", s = 60)
    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(7.2), color = "grey", label ="$\\beta=7.20$", s = 60)

    x_low, x_high = [1e-3,2]
    plt.xlim([x_low,x_high])
    plt.ylim([-4,0.1])


    plt.xscale("log")
    plt.legend()
    plt.xlabel("$\\frac{P^2}{m_\pi^2}$")
    plt.ylabel("$\\frac{P \cot(\delta)}{m_\pi}$")
    plt.grid()
    plt.savefig("plots/P_cot_all.pdf", bbox_inches="tight")
    xarr = np.logspace(np.log(x_low), np.log(x_high), 100)
    for i in range(len(a2)):
        plt.plot(xarr, fit_func_p2(xarr, a2[i], b2[i]), color = "maroon", lw = 0.5, alpha = 0.8)#, ls = "dashed")
    plt.savefig("plots/P_cot_all_fit.pdf", bbox_inches="tight")
    plt.show()

def plot_tan_PS():
    plt.figure(figsize=(8,4.8))
    P_cot_PS = []
    P_cot_PS_p = []
    P_cot_PS_m = []
    P_2 = []
    P_2_p = []
    P_2_m = []
    a2 = []
    b2 = []
    m_pi_rho_arr = []
    beta_arr = []
    m_arr = []
    all_meas = errcl.measurement("all_SP4")
    all_meas.read_from_HDF()

    for result_name in all_meas.result_names:
        search_str = "phase_shift_Fabian_level_0_Scattering_I2_SP(4)"
        # print(result_name)
        if result_name[:len(search_str)] == search_str:
            search_str = "tan_PS/P_prime"
            if result_name[len(result_name)-len(search_str):] == search_str:
                info_str = result_name[:len(result_name)-len(search_str)]
                P_cot_PS.append(all_meas.results[result_name].median)
                P_cot_PS_p.append(all_meas.results[result_name].ep)
                P_cot_PS_m.append(all_meas.results[result_name].em)
                P_2.append(all_meas.results[result_name[:len(result_name)-len(search_str)]+"P_2_prime"].median)
                P_2_p.append(all_meas.results[result_name[:len(result_name)-len(search_str)]+"P_2_prime"].ep)
                P_2_m.append(all_meas.results[result_name[:len(result_name)-len(search_str)]+"P_2_prime"].em)
                beta = all_meas.infos[info_str+"beta"]
                beta_arr.append(beta)
                m = all_meas.infos[info_str+"m_1"]
                m_arr.append(m)
                m_pi_rho = all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].median[0]
                m_pi_rho_arr.append(m_pi_rho)
                # print("beta = %1.3f, m = %1.3f, m_pi_rho = %1.3f"%(beta, m, m_pi_rho))
        search_str = "phase_shift_fit_P_cot_PS_SP(4)"
        if result_name[:len(search_str)] == search_str:
            search_str = "a2"
            if result_name[len(result_name)-len(search_str):] == search_str:
                info_str = result_name[:len(result_name)-len(search_str)]
                a2.append(all_meas.results[result_name].median[0])
                b2.append(all_meas.results[result_name[:len(result_name)-len(search_str)]+"b2"].median[0])

            
    norm = matplotlib.colors.Normalize(vmin=min(m_pi_rho_arr), vmax=max(m_pi_rho_arr), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    time_color = np.array([(mapper.to_rgba(v)) for v in m_pi_rho_arr])
    plt.colorbar(mappable=mapper, label = "$\\frac{m_\pi}{m_\\rho}$")
    beta_used = []
    for i in range(len(P_2)):   
        plt.errorbar(x=P_2[i], xerr=[P_2_m[i],P_2_p[i]], y=P_cot_PS[i], yerr=[P_cot_PS_m[i],P_cot_PS_p[i]], marker = marker_beta(beta_arr[i]), color = time_color[i], capsize=5, markersize=10)

    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(6.9), color = "grey", label ="$\\beta=6.90$", s = 60)
    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(7.05), color = "grey", label ="$\\beta=7.05$", s = 60)
    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(7.2), color = "grey", label ="$\\beta=7.20$", s = 60)

    x_low, x_high = [1e-3,0.7]
    plt.xlim([x_low,x_high])
    plt.ylim([-1.7,0])
    plt.xscale("log")
    plt.legend()
    plt.xlabel("$\\frac{P^2}{m_\pi^2}$")
    plt.ylabel("$\\frac{\\tan(\delta_0)m_\pi}{P}$")
    plt.grid()
    # plt.grid(ls = "dashed", color = "black")
    plt.savefig("plots/tan_PS_all.pdf", bbox_inches="tight")
    xarr = np.logspace(np.log(x_low), np.log(1), 1000)
    for i in range(len(a2)):
        # print(a2[i],b2[i])
        yarr = [fit_func_tan_p2(x, a2[i], b2[i]) for x in xarr]
        # print(xarr)
        # print(yarr)
        plt.plot(xarr, yarr, color = "maroon", lw = 0.5, alpha = 0.8)#, ls = "dashed")
    plt.savefig("plots/tan_PS_all_fit.pdf", bbox_inches="tight")
    plt.show()

def plot_PS():
    P_cot_PS = []
    P_cot_PS_p = []
    P_cot_PS_m = []
    P_2 = []
    P_2_p = []
    P_2_m = []
    P = []
    a2 = []
    b2 = []
    m_pi_rho_arr = []
    beta_arr = []
    m_arr = []
    all_meas = errcl.measurement("all_SP4")
    all_meas.read_from_HDF()

    for result_name in all_meas.result_names:
        search_str_1 = "phase_shift_Fabian_level_0_Scattering_I2_SP(4)"
        # print(result_name)
        if result_name[:len(search_str_1)] == search_str_1:
            search_str_2 = "/PS"
            if result_name[len(result_name)-len(search_str_2):] == search_str_2:
                info_str = result_name[:len(result_name)-len(search_str_2)]
                P_cot_PS.append(all_meas.results[result_name].median)
                P_cot_PS_p.append(all_meas.results[result_name].ep)
                P_cot_PS_m.append(all_meas.results[result_name].em)
                P_2.append(all_meas.results[result_name[:len(result_name)-len(search_str_2)]+"/s_pipi_prime"].median)
                P_2_p.append(all_meas.results[result_name[:len(result_name)-len(search_str_2)]+"/s_pipi_prime"].ep)
                P_2_m.append(all_meas.results[result_name[:len(result_name)-len(search_str_2)]+"/s_pipi_prime"].em)
                beta = all_meas.infos[info_str+"/beta"]
                beta_arr.append(beta)
                m = all_meas.infos[info_str+"/m_1"]
                m_arr.append(m)
                m_pi_rho = all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].median[0]
                m_pi_rho_arr.append(m_pi_rho)
                # print("beta = %1.3f, m = %1.3f, m_pi_rho = %1.3f"%(beta, m, m_pi_rho))
                search_str_3 = "phase_shift_fit_P_cot_PS_SP(4)"
                # print("phase_shift_fit_P_cot_PS_SP(4)"+result_name[len(search_str_1):len(result_name)-len(search_str_2)-5]+"/a2")
                # exit()
                # if ("phase_shift_fit_P_cot_PS_SP(4)"+result_name[:len(result_name)-search_str_3-5]+"/a2") in all_meas.result_names:
                #     search_str = "a2"
                    # if result_name[len(result_name)-len(search_str):] == search_str:
                    #     info_str = result_name[:len(result_name)-len(search_str)]
                    #     a2.append(all_meas.results[result_name].median[0])
                    #     b2.append(all_meas.results[result_name[:len(result_name)-len(search_str)]+"b2"].median[0])

            
    norm = matplotlib.colors.Normalize(vmin=min(m_pi_rho_arr), vmax=max(m_pi_rho_arr), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    time_color = np.array([(mapper.to_rgba(v)) for v in m_pi_rho_arr])
    plt.colorbar(mappable=mapper, label = "$\\frac{m_\pi}{m_\\rho}$")
    beta_used = []
    for i in range(len(P_2)):   
        plt.errorbar(x=P_2[i], xerr=[P_2_m[i],P_2_p[i]], y=P_cot_PS[i], yerr=[P_cot_PS_m[i],P_cot_PS_p[i]], marker = marker_beta(beta_arr[i]), color = time_color[i], capsize=5, markersize=10)

    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(6.9), color = "grey", label ="$\\beta=6.90$", s = 60)
    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(7.05), color = "grey", label ="$\\beta=7.05$", s = 60)
    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(7.2), color = "grey", label ="$\\beta=7.20$", s = 60)

    plt.xlim([4,6.2])
    plt.ylim([-1,0])
    plt.xscale("log")
    plt.legend()
    plt.xlabel("$\\frac{s}{m_\\pi^2}$")
    plt.ylabel("$\delta$")
    plt.grid()
    plt.savefig("plots/PS.pdf", bbox_inches="tight")
    xarr = np.logspace(np.log(1.1), np.log(10), 100)
    for i in range(len(a2)):
        yarr = [fit_func_PS_s(x, a2[i], b2[i]) for x in xarr]
        plt.plot(xarr, yarr, color = "grey", lw = 0.5)
    plt.savefig("plots/PS_fit.pdf", bbox_inches="tight")
    plt.show()

def trivial_energy_level(N_L_inv, level, m):
    def E_pipi(m, P):
        return 2*np.arccosh(np.cosh(m)+2*np.sin(P/2)**2)
        # return 2*np.sqrt(m**2+P**2)
    return E_pipi(m, level*2*np.pi*N_L_inv)

def plot_energy_levels():
    E_pipi = []
    E_pipi_p = []
    E_pipi_m = []
    m_inf = []
    m_inf_p = []
    m_inf_m = []
    N_L_inv = []
    beta_m_arr = []
    beta_arr = []
    m_arr = []

    all_meas = errcl.measurement("all_SP4")
    all_meas.read_from_HDF()
    for result_name in all_meas.result_names:
        search_str = "infinite_volume_Fabian"
        if result_name[:len(search_str)] == search_str:
            search_str = "_pipi/E"
            if result_name[len(result_name)-len(search_str):] == search_str:
                beta_m = result_name[55:60]+result_name[63:69]
                level = int(result_name[29])
                print(result_name, level)
                if not (beta_m in beta_m_arr):
                    beta_m_arr.append(beta_m)
                    print(beta_m)
                    beta_arr.append(float(beta_m[:5]))
                    m_arr.append(float(beta_m[5:]))
                    E_pipi.append([None, None])
                    E_pipi_p.append([None, None])
                    E_pipi_m.append([None, None])
                    N_L_inv.append([]) 
                    m_inf.append([]) 
                    m_inf_p.append([]) 
                    m_inf_m.append([]) 
                E_pipi[beta_m_arr.index(beta_m)][level] = all_meas.results[result_name].median
                E_pipi_p[beta_m_arr.index(beta_m)][level] = all_meas.results[result_name].ep
                E_pipi_m[beta_m_arr.index(beta_m)][level] = all_meas.results[result_name].em
                N_L_inv[beta_m_arr.index(beta_m)] = all_meas.results[result_name[:len(result_name)-2]+"/N_Ls_inv"].median
                m_inf[beta_m_arr.index(beta_m)] = all_meas.results[result_name[:len(result_name)-6]+"pi/m_inf"].median[0]
                m_inf_p[beta_m_arr.index(beta_m)] = all_meas.results[result_name[:len(result_name)-6]+"pi/m_inf"].median_p[0]
                m_inf_m[beta_m_arr.index(beta_m)] = all_meas.results[result_name[:len(result_name)-6]+"pi/m_inf"].median_m[0]
                    

    for i in range(len(E_pipi)):
        xarr = np.linspace(0,1/7., 500)
        yarr1 =  trivial_energy_level(xarr, 1, m_inf[i])
        yarr2 =  trivial_energy_level(xarr, 2, m_inf[i])
        print(N_L_inv[i], E_pipi[i])
        plt.errorbar(x=N_L_inv[i], y=E_pipi[i][0],yerr=[E_pipi_p[i][0],E_pipi_m[i][0]], ls = "", capsize=5, markersize=10, marker = "^", color = "purple", label = "$E_{\pi\pi}$ (n=0)")
        if E_pipi[i][1] != None:
            plt.errorbar(x=N_L_inv[i], y=E_pipi[i][1],yerr=[E_pipi_p[i][1],E_pipi_m[i][1]], ls = "", capsize=5, markersize=10, marker = "v", color = "purple", alpha = 0.5, label = "$E_{\pi\pi}$ (n=1)")
        plt.axhline(2*m_inf[i], color = "blue", label = "$2m_\pi$")
        plt.axhline(4*m_inf[i], color = "green", label = "$4m_\pi$")
        plt.fill_between(x = [0,1/7.], y1 = 2*m_inf_m[i], y2 = 2*m_inf_p[i], color = "blue", alpha = 0.5)
        plt.fill_between(x = [0,1/7.], y1 = 4*m_inf_m[i], y2 = 4*m_inf_p[i], color = "green", alpha = 0.5)
        plt.title(beta_m_arr[i])
        plt.xlabel("$N_L^{-1}$")
        plt.ylabel("E")
        plt.plot(xarr, yarr1, color = "black", ls="dotted", label="level 1")
        plt.plot(xarr, yarr2, color = "black", ls="--", label="level 2")
        plt.legend()
        plt.grid()
        plt.savefig("plots/trivial_%s.pdf"%beta_m_arr[i])
        # plt.show()
        plt.clf()


def sigma_of_P(P, coeffs, mass):                                              # P - momentum, coeffs contains coeffs of expansions
    E = 2*np.arccosh(np.cosh(mass)+2*np.sin(P/2)**2)
    f = 0
    for i in range(len(coeffs)):
        if i == 0:
            f += coeffs[0]
        else:
            f += coeffs[i]*P**(2*i)
    
    return 16*np.pi/(E**2*(1+(f/P)**2))

def plot_sigma():
    sigma = []
    sigma_p = []
    sigma_m = []
    P = np.linspace(0,2,50)
    # v = []
    # E = []
    # s = []
    m_pi_rho_arr = []
    beta_arr = []
    m_arr = []
    all_meas = errcl.measurement("all_SP4")
    all_meas.read_from_HDF()
    print("all meas read")

    ind = 0

    for result_name in all_meas.result_names:
        search_str = "phase_shift_fit_P_cot_PS_SP(4)"
        if result_name[:len(search_str)] == search_str:
            search_str = "a2"
            sigma.append([])
            sigma_p.append([])
            sigma_m.append([])
            if result_name[len(result_name)-len(search_str):] == search_str:
                search_str = "lim_0/"                                               # which limit for P do you want?
                if search_str in result_name:
                    info_str = result_name[:len(result_name)-2]
                    beta = all_meas.infos[info_str+"beta"]
                    beta_arr.append(beta)
                    m = all_meas.infos[info_str+"m_1"]
                    m_arr.append(m)
                    m_pi_rho = all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].median[0]
                    m_pi_rho_arr.append(m_pi_rho)

                    mass = all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/mass_Goldstone"%(beta,m,m)].median[0]
                    coeffs = []
                    coeffs.append(all_meas.results[result_name[:len(result_name)-2]+"a2"].median[0])
                    coeffs.append(all_meas.results[result_name[:len(result_name)-2]+"b2"].median[0])
                    for i in range(len(P)):
                        sigma[ind].append(sigma_of_P(P[i], coeffs, mass))
                    plt.plot(P, sigma[ind])

                    ind += 1
                    print(result_name)
    plt.show()
    #             else:
    #                 exit()



    #             info_str = result_name[:len(result_name)-len(search_str)]
    #             P_cot_PS.append(all_meas.results[result_name].median)
    #             P_cot_PS_p.append(all_meas.results[result_name].ep)
    #             P_cot_PS_m.append(all_meas.results[result_name].em)
    #             P_2.append(all_meas.results[result_name[:len(result_name)-len(search_str)]+"P_2_prime"].median)
    #             P_2_p.append(all_meas.results[result_name[:len(result_name)-len(search_str)]+"P_2_prime"].ep)
    #             P_2_m.append(all_meas.results[result_name[:len(result_name)-len(search_str)]+"P_2_prime"].em)
    #             beta = all_meas.infos[info_str+"beta"]
    #             beta_arr.append(beta)
    #             m = all_meas.infos[info_str+"m_1"]
    #             m_arr.append(m)
    #             m_pi_rho = all_meas.results["infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi/m_pi_rho"%(beta,m,m)].median[0]
    #             m_pi_rho_arr.append(m_pi_rho)
    #             # print("beta = %1.3f, m = %1.3f, m_pi_rho = %1.3f"%(beta, m, m_pi_rho))
    #     search_str = "phase_shift_fit_P_cot_PS_SP(4)"
    #     if result_name[:len(search_str)] == search_str:
    #         search_str = "a2"
    #         if result_name[len(result_name)-len(search_str):] == search_str:
    #             info_str = result_name[:len(result_name)-len(search_str)]
    #             a2.append(all_meas.results[result_name].median)
    #             b2.append(all_meas.results[result_name[:len(result_name)-len(search_str)]+"b2"].median)


    # norm = matplotlib.colors.Normalize(vmin=min(m_pi_rho_arr), vmax=max(m_pi_rho_arr), clip=True)
    # mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    # time_color = np.array([(mapper.to_rgba(v)) for v in m_pi_rho_arr])
    # plt.colorbar(mappable=mapper, label = "$\\frac{m_\pi}{m_\\rho}$")
    # beta_used = []
    # for i in range(len(P_2)):   
    #     plt.errorbar(x=P_2[i], xerr=[P_2_m[i],P_2_p[i]], y=P_cot_PS[i], yerr=[P_cot_PS_m[i],P_cot_PS_p[i]], marker = marker_beta(beta_arr[i]), color = time_color[i], capsize=5, markersize=10)


    # plt.scatter((-10, -9), y = (0,0), marker = marker_beta(6.9), color = "grey", label ="$\\beta=6.90$", s = 60)
    # plt.scatter((-10, -9), y = (0,0), marker = marker_beta(7.05), color = "grey", label ="$\\beta=7.05$", s = 60)
    # plt.scatter((-10, -9), y = (0,0), marker = marker_beta(7.2), color = "grey", label ="$\\beta=7.20$", s = 60)

    # x_low, x_high = [1e-3,2]
    # plt.xlim([x_low,x_high])
    # plt.ylim([-4,0.1])


    # plt.xscale("log")
    # plt.legend()
    # plt.xlabel("$\\frac{P^2}{m_\pi^2}$")
    # plt.ylabel("$\\frac{P \cot(\delta)}{m_\pi}$")
    # plt.grid()
    # plt.savefig("plots/P_cot_all.pdf", bbox_inches="tight")
    # xarr = np.logspace(np.log(x_low), np.log(x_high), 100)
    # for i in range(len(a2)):
    #     plt.plot(xarr, fit_func_p2(xarr, a2[i], b2[i]), color = "maroon", lw = 0.5, alpha = 0.8)#, ls = "dashed")
    # plt.savefig("plots/P_cot_all_fit.pdf", bbox_inches="tight")
    # plt.show()




def main():

    #######################################################

    # add_m_pi_m_rho()
    # plot_a_0_vs_m_pi_m_rho()
    # plot_tan_PS()
    # plot_P_cot_PS()
    # plot_PS()
    # plot_energy_levels()
    # plot_a_0_vs_m_f_pi(squared = False)
    # plot_a_0_vs_m_f_pi(squared = True)

    plot_sigma()

    #######################################################


if __name__ == "__main__":
    main()



