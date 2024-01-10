import os
import matplotlib.pyplot as plt
import read_HDF5_logfile as HDF_log
import error_classes as errcl 
import numpy as np
import matplotlib
from matplotlib import cm
import math

color_arr = ["blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick"]

font = {'size'   : 16}
matplotlib.rc('font', **font)

def plot(filelist, x, y, error = False, fit = False, zoom = None, title = "", save_app = "", logx = False, logy = False, color = 100*["red",], save = True, show = True):
    zoom_txt = ""
    log_txt = ""
    fit_txt = ""

    if error:
        plot_errorbar(filelist, x, y, color)
    else:
        plot_scatter(filelist, x, y, color)
    if fit:
        fit_txt = "_fit"
        plot_fit(filelist, x, y, error, color[0])

    log_axis(logx, logy)

    if save_app != "":
        save_app = "_"+save_app
    if zoom != None:
        zoom_txt = "_zoom"
        plt.xlim((zoom[0], zoom[1]))
        plt.ylim((zoom[2], zoom[3]))
    plt.xlabel(x)
    plt.ylabel(y)
    plt.legend()
    plt.title(title)
    plt.grid()

    y = y.replace("/", "_")

    if save:
        plt.savefig("plots/"+y+"_vs_"+x+log_txt+save_app+zoom_txt+".pdf")
    if show:
        plt.show()

def log_axis(logx, logy):
    if logx and logy:
        plt.xscale("log")
        plt.yscale("log")
        log_txt = "_logxy"
    if logx:
        plt.xscale("log")
        log_txt = "_logx"
    if logy:
        plt.yscale("log")
        log_txt = "_logy"

def plot_scatter(filelist, x, y, color = 100*["red",]):
    for i in range(len(filelist)):
        meas = errcl.measurement(filelist[i])
        meas.read_from_HDF()
        plt.scatter(x=meas.results[x].median, y=meas.results[y].median, color = color[i]) 
        
def plot_errorbar(filelist, x, y, color = 100*["red",]):
    for i in range(len(filelist)):
        meas = errcl.measurement(filelist[i])
        meas.read_from_HDF()
        # if color[i] == "red":
        #     if meas.results[y].median[0] > 0:
        #         color[i] = "orange"
        #     if meas.results[y].median_m[0] > 0 and meas.results[y].median_p[0] > 0:
        #         color[i] = "magenta"
        # if color[i] == "blue":
        #     if meas.results[y].median[0] > 0:
        #         color[i] = "c"
        #     if meas.results[y].median_m[0] > 0 and meas.results[y].median_p[0] > 0:
        #         color[i] = "green"
        plt.errorbar(x=meas.results[x].median,xerr=[meas.results[x].em,meas.results[x].ep], y=meas.results[y].median[0],yerr=[meas.results[y].em,meas.results[y].ep], color = color[i], ls = "", capsize=5, markersize=10) 

def plot_fit(filelist, x, y, error, color = "red"):
    if x == "P_2_prime":
        if y == "tan_PS/P_prime":
            fit_tan_PS(filelist, x, y, error, color)
        elif y == "P_cot_PS_prime":
            fit_P_cot_PS(filelist, x, y, error, color)

def fit_func_p2(P2, a, b):
    return a + b*P2
def fit_func_tan_p2(P2, a, b):
    res = 1/(a + b*P2)
    if res > 0:
        return np.nan
    return res
def fit_func_PS_s(s, a, b):
    P2 = s-1
    return np.arctan(1/(a + b*P2))
# def fit_func_p4(P2, a, b, c):
#     return a + b*P2 + c*P2*P2

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
        

# def fit_func_tan_P2(P2, a, b):
#     return 1/(b*P2+a)
# def fit_func_tan_P4(P2, a, b, c):
#     return 1/(a+b*P2+c*P2*P2)

# def fit_tan_PS(filelist, x, y, error, color = "red"):
#     bm_str = filelist[0][41:65]
#     meas = errcl.measurement("phase_shift_fit_P_cot_PS_"+bm_str)
#     if meas.file_exists():
#         meas.read_from_HDF()
#         fit_result = meas.results
#         num = len(meas.results[x].median)
#         # xarr = np.linspace(min(meas.results[x].median)*0.9,max(meas.results[x].median)*1.1, 500)
#         xarr = np.linspace(0,2, 500)
#         yarr_2 = fit_func_tan_P2(xarr, fit_result["a2"].median[0], fit_result["b2"].median[0])
#         plt.plot(xarr, yarr_2, color = color, label = "O(P2)")
#         if error:
#             yarr_2m = fit_func_tan_P2(xarr, fit_result["a2"].median_m[0], fit_result["b2"].median[0])  
#             yarr_2p = fit_func_tan_P2(xarr, fit_result["a2"].median_p[0], fit_result["b2"].median[0])  
#             plt.fill_between(x=xarr, y1=yarr_2m, y2=yarr_2p, color = color, alpha = 0.3)
#         if num > 2:
#             yarr_4 = fit_func_tan_P4(xarr, fit_result["a4"].median[0], fit_result["b4"].median[0], fit_result["c4"].median[0])
#             plt.plot(xarr, yarr_4, color = color, ls = "--", label = "O(P4)")
#             if error:
#                 yarr_4m = fit_func_tan_P4(xarr, fit_result["a4"].median_m[0], fit_result["b4"].median[0], fit_result["c4"].median[0])  
#                 yarr_4p = fit_func_tan_P4(xarr, fit_result["a4"].median_p[0], fit_result["b4"].median[0], fit_result["c4"].median[0])  
#                 plt.fill_between(x=xarr, y1=yarr_4m, y2=yarr_4p, color = color, alpha = 0.1)



    # xarr = np.logspace(np.log(min(x)[0]*0.9),np.log(max(x)[0]*1.1), 500)
    # yarr_2 = fit_func_p2(xarr, fit_result["a2"].median[0], fit_result["b2"].median[0])
    # yarr_2p = fit_func_p2(xarr, fit_result["a2"].median_p[0], fit_result["b2"].median[0])
    # yarr_2m = fit_func_p2(xarr, fit_result["a2"].median_m[0], fit_result["b2"].median[0])


    # plt.plot(xarr, yarr_4, color = color_arr[1], label = "$\\mathcal{O}(P^4)$") #: a_0 = %1.3f +  %1.3f -  %1.3f$"%(fit_result["a_0_4"].median[0],fit_result["a_0_4"].ep[0],fit_result["a_0_4"].em[0]))
    # plt.fill_between(x=xarr, y1=yarr_4m, y2=yarr_4p, color = color_arr[1], alpha = 0.3)

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
    plt.ylim([-1.02,0])
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





def main():

    #######################################################

    # add_m_pi_m_rho()
    # plot_a_0_vs_m_pi_m_rho()
    plot_tan_PS()
    plot_P_cot_PS()
    # plot_PS()
    # plot_energy_levels()
    # plot_a_0_vs_m_f_pi(squared = False)
    # plot_a_0_vs_m_f_pi(squared = True)

    #######################################################




    # filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot_SP4_lev0", "str")
    # print(filelist)
    # plot(filelist, "P_2_prime", "P_cot_PS_prime", error=True, log=True)
    # plot(filelist, "P_2_prime", "tan_PS/P_prime", error=True, log=True)

    #######################################################
    
    # plt.fill_between((0.01,2), (-1, -1), (-0.095900719, -0.095900719), color ="blue", alpha = 0.3)
    # plt.fill_between((0.01,2), (-0.095900719, -0.095900719),(0, 0), color ="red", alpha = 0.3)
    # plt.fill_between((0.01,2), (0,0),(0.472894247, 0.472894247), color ="blue", alpha = 0.3)
    # plt.fill_between((0.01,2), (0.472894247,0.472894247),(1,1), color ="red", alpha = 0.3)
    # plt.fill_between((0.01,2), (1,1),(1.441591313,1.441591313), color ="blue", alpha = 0.3)
    # plt.fill_between((0.01,2), (1.441591313,1.441591313),(2,2), color ="red", alpha = 0.3)
    # plot(filelist, "P_prime", "q_2", error = True, logx = True, zoom=[0.01, 2, -0.5, 1.8], color = 100*["black",], title = "red => Zeta(q^2) = positive, blue => Zeta(q^2) = negative")
    
    #######################################################

    #######################################################

    # clr = 100*["red",]
    # i = 0
    # with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot_SP4_lev01_split", "r") as f:
    #     for lines in f.readlines():
    #         filelist = lines.split()
    #         bm_str = filelist[0][41:65]
    #         if int(filelist[0][25]) == 1:
    #             clr = 100*["blue",]
    #         else: 
    #             clr = 100*["red",]
    #         plot(filelist, "P_2_prime", "P_cot_PS_prime",show = False, save = False, fit = False, error=True, logx=True, save_app=bm_str, zoom=[1e-4, 2, -5,5], color = clr)
 
    #         i += 1
    # # plt.grid()
    # plt.savefig("plots/all_P_cot_PS_lv1.pdf")
    # plt.show()

    # i = 0
    # with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot_SP4_lev0_split", "r") as f:
    #     for lines in f.readlines():
    #         filelist = lines.split()
    #         bm_str = filelist[0][41:65]
    #         plot(filelist, "P_2_prime", "P_cot_PS_prime",show = False, save = False, fit = False, error=True, logx=True, save_app=bm_str, color=100*[color_arr[i],], zoom=[1e-4, 2, -5,5])
 
    #         i += 1
    # plt.grid()
    # plt.savefig("plots/all_P_cot_PS.pdf")
    # plt.show()

    # i = 0
    # with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot_SP4_lev0_split", "r") as f:
    #     for lines in f.readlines():
    #         filelist = lines.split()
    #         bm_str = filelist[0][41:65]
    #         plot(filelist, "P_2_prime", "tan_PS/P_prime",show = False, save = False, fit = False, error=True, logx=True, title = bm_str, save_app=bm_str, color=100*[color_arr[i],], zoom=[1e-4, 2, -5,5])
    #         i += 1
    # plt.grid()
    # plt.savefig("plots/all_tan_PS.pdf")
    # plt.show()

    #################################

    # filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_fit_all", "str")
    # plot(filelist, "P_2_prime", "P_cot_PS_prime", error=True, log=True)
    


if __name__ == "__main__":
    main()



    
    # results = get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot")
    # # results = get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot_SU3")
    # # fit_result = get_result_file("phase_shift_fit_P_cot_PS")

    # P_cot_PS, P_cot_PS_err, P_2, P_2_err, rho_pi = [[],[],[],[],[]]
    # for i in range(len(results)):
    #     rho_pi.append(results[i]["infos/m_rho_pi"])
    #     P_2.append(results[i]["P_2"].median)
    #     print("E_pipi: ", results[i]["E_pipi"].median)
    #     print("tan_PS: ", results[i]["tan_PS"].median)
    #     print("cot_PS: ", results[i]["cot_PS"].median)
    #     print("P: ", results[i]["P"].median)
    #     print("P_2: ", results[i]["P_2"].median)
    #     print("m_inf: ", results[i]["mass_Goldstone"].median)
    #     print("N_L: ", results[i]["N_L"].median)
    #     print("\n")
    #     P_2_err.append([results[i]["P_2"].ep,results[i]["P_2"].em])
    #     P_cot_PS.append(results[i]["cot_PS"].median[0])
    #     P_cot_PS_err.append([results[i]["cot_PS"].ep,results[i]["cot_PS"].em])
    
    # norm = matplotlib.colors.Normalize(vmin=min(rho_pi), vmax=max(rho_pi), clip=True)
    # mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    # time_color = np.array([(mapper.to_rgba(v)) for v in rho_pi])
    # plt.colorbar(mappable=mapper, label = "$\\frac{m_{\\rho}}{m_{\\pi}}$")
    # for i in range(len(results)):
    #     # plt.errorbar(x=P_2[i], xerr = P_2_err[i], y=P_cot_PS[i], yerr=P_cot_PS_err[i], color = time_color[i], ls = "", capsize=5, markersize=10)
    #     plt.scatter(x=P_2[i], y=P_cot_PS[i], color = time_color[i])
    # if fit:
    #     plot_fits_P_cot_PS(P_2, fit_result)


    # # plt.axvline(0.4, color = "black", ls="--")
    # log_txt = ""
    # if log:
    #     plt.xscale("log")
    #     log_txt = "_log"
    # plt.ylabel("$\\frac{Pcot(\delta)}{m_{\pi}}$")
    # plt.xlabel("$\\frac{P^2}{m_{\pi}^2}$")
    # plt.legend(fontsize = "x-small")
    # plt.title(title)
    # plt.grid()
    # zoom_txt = ""
    # if zoom == plot_tan_PS.__defaults__[0]:
    #     zoom_txt = "_zoom"
    # if zoom != False:
    #     plt.xlim((zoom[0], zoom[1]))
    #     plt.ylim((zoom[2], zoom[3]))
    # plt.savefig("plots/P_cot_PS"+log_txt+save_app+zoom_txt+".pdf")
    # plt.show()
    # plt.clf()


























































color_arr = ["blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick","blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick","blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick",]

def create_all_filenames():
    PATH = "output/result_files/"
    temp = "phase_shift_Fabian"
    filelist = os.listdir(PATH)
    resultfile_list = []
    num = len(temp)
    for file in filelist:
        length = len(file)
        if file[:num] == temp:
            resultfile_list.append(file[:length-5])         

    with open("input/filenames_phase_shift_plot_all", "w") as file:
        for filename in resultfile_list:
            file.write("%s"%filename+"\n")

def get_results_from_files(filelist):
    filelist = np.genfromtxt(filelist, "str")
    results = []
    for i in range(len(filelist)):
        print(filelist[i])
        results.append(get_result_file(filelist[i]))
    return results

def get_result_file(file):
    print(file)
    dic = {}
    meas = errcl.measurement(file)
    meas.read_from_HDF()
    meas.print_everything()
    for key in meas.result_names:
        # print(key)
        dic[key] = meas.results[key]
    for key in meas.infos.keys():
        dic["infos/"+key] = meas.infos[key]
    return dic


def plot_q_P():
    import generalizedzeta as gz
    results = get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot")
    
    P, P_err, q, q_err, rho_pi = [[],[],[],[],[]]
    for i in range(len(results)):
        rho_pi.append(results[i]["infos/m_rho_pi"])
        P.append(results[i]["P"].median)
        P_err.append([results[i]["P"].ep,results[i]["P"].em])
        q.append(results[i]["q"].median[0])
        q_err.append([results[i]["q"].ep,results[i]["q"].em])
    # xarr = np.linspace(min(q_2)*0.9, max(q_2)*1.1, 500)
    # xarr = np.linspace(-1, 1.5, 500)
    # yarr = []
    # for x in xarr:
    #     yarr.append(gz.Zeta(x))
    # plt.plot(xarr, yarr)
    plt.ylabel("P")
    plt.xlabel("q")

    norm = matplotlib.colors.Normalize(vmin=min(rho_pi), vmax=max(rho_pi), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    time_color = np.array([(mapper.to_rgba(v)) for v in rho_pi])
    plt.colorbar(mappable=mapper, label = "$\\frac{m_{\\rho}}{m_{\\pi}}$")
    for i in range(len(results)):
        plt.errorbar(x=q[i], xerr = q_err[i], y=P[i], yerr=P_err[i], color = time_color[i], ls = "", capsize=5, markersize=10)
    plt.grid()
    # plt.xlim((-0.5,1.5))
    # plt.ylim((-55,10))
    plt.savefig("q_vs_P.pdf")
    plt.show()

def plot_dispersion_relation():
    import phase_shift
    results = get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot")
    
    P_2_prime, P_2_prime_err, s_pipi_prime, s_pipi_prime_err, mass_Goldstone, rho_pi = [[],[],[],[],[],[]]
    P, P_err, E_pipi, E_pipi_err = [[],[],[],[]]
    for i in range(len(results)):
        print(results[i]["E_pipi_prime"].median[0])
        if results[i]["infos/beta"] == 6.9 and results[i]["infos/m_1"] == -0.90:
            mass_Goldstone.append(results[i]["mass_Goldstone"].median)
            rho_pi.append(results[i]["infos/m_rho_pi"])
            P.append(results[i]["P"].median)
            P_err.append([results[i]["P"].ep,results[i]["P"].em])
            E_pipi.append(results[i]["E_pipi"].median[0])
            E_pipi_err.append([results[i]["E_pipi"].ep,results[i]["E_pipi"].em])
            P_2_prime.append(results[i]["P_2_prime"].median)
            P_2_prime_err.append([results[i]["P_2_prime"].ep,results[i]["P_2_prime"].em])
            s_pipi_prime.append(results[i]["s_pipi_prime"].median[0])
            s_pipi_prime_err.append([results[i]["s_pipi_prime"].ep,results[i]["s_pipi_prime"].em])

    xarr = np.linspace(min(E_pipi)*0.9, max(E_pipi)*1.1, 500)
    # xarr2 = []
    yarr_lat = []
    yarr_classic = []
    for x in xarr:
        # xarr2.append(x*x)
        yarr_lat.append(phase_shift.generalized_mom(x, np.mean(mass_Goldstone)))
        yarr_classic.append(phase_shift.generalized_mom_non_lattice(x, np.mean(mass_Goldstone)))
    plt.plot(xarr, yarr_lat, label="disp rel lattice")
    plt.plot(xarr, yarr_classic, label="disp rel classic")

    norm = matplotlib.colors.Normalize(vmin=min(rho_pi), vmax=max(rho_pi), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    time_color = np.array([(mapper.to_rgba(v)) for v in rho_pi])
    plt.colorbar(mappable=mapper, label = "$\\frac{m_{\\rho}}{m_{\\pi}}$")
    for i in range(len(E_pipi)):
        plt.errorbar(x=E_pipi[i], xerr = E_pipi_err[i], y=P[i], yerr=P_err[i], color = time_color[i], ls = "", capsize=5, markersize=10)
    plt.grid()
    # plt.xlim((-0.5,1.5))
    # plt.ylim((-55,10))
    plt.legend()
    plt.savefig("dispersion_relation.pdf")
    plt.show()



# def fit_func_p2(P2, a, b):
#     return a + b*P2
# def fit_func_p4(P2, a, b, c):
#     return a + b*P2 + c*P2*P2
# def fit_func_p6(P2, a, b, c, d):
#     return a + b*P2 + c*P2*P2 + d*P2*P2*P2
# def fit_func_p2_fixed_400(P2, b):
#     return -1/400 + b*P2
# def fit_func_p4_fixed_400(P2, b, c):
#     return -1/400 + b*P2 + c*P2*P2

# def plot_fits_P_cot_PS(P_2, fit_result):
#     xarr = np.logspace(np.log(min(P_2)[0]*0.9),np.log(max(P_2)[0]*1.1), 500)
#     yarr_2 = fit_func_p2(xarr, fit_result["a2"].median[0], fit_result["b2"].median[0])
#     yarr_2p = fit_func_p2(xarr, fit_result["a2"].median_p[0], fit_result["b2"].median[0])
#     yarr_2m = fit_func_p2(xarr, fit_result["a2"].median_m[0], fit_result["b2"].median[0])
#     yarr_2_fixed = fit_func_p2_fixed_400(xarr, fit_result["b2_fixed"].median[0])
#     yarr_2p_fixed = fit_func_p2_fixed_400(xarr, fit_result["b2_fixed"].median[0])
#     yarr_2m_fixed = fit_func_p2_fixed_400(xarr, fit_result["b2_fixed"].median[0])
#     yarr_4 = fit_func_p4(xarr, fit_result["a4"].median[0], fit_result["b4"].median[0], fit_result["c4"].median[0])
#     yarr_4p = fit_func_p4(xarr, fit_result["a4"].median_p[0], fit_result["b4"].median[0], fit_result["c4"].median[0])
#     yarr_4m = fit_func_p4(xarr, fit_result["a4"].median_m[0], fit_result["b4"].median[0], fit_result["c4"].median[0])
#     yarr_4_fixed = fit_func_p4_fixed_400(xarr, fit_result["b4_fixed"].median[0], fit_result["c4_fixed"].median[0])
#     yarr_4p_fixed = fit_func_p4_fixed_400(xarr, fit_result["b4_fixed"].median_p[0], fit_result["c4_fixed"].median[0])
#     yarr_4m_fixed = fit_func_p4_fixed_400(xarr, fit_result["b4_fixed"].median_m[0], fit_result["c4_fixed"].median[0])
#     yarr_6 = fit_func_p6(xarr, fit_result["a6"].median[0], fit_result["b6"].median[0], fit_result["c6"].median[0], fit_result["d6"].median[0])
#     yarr_6p = fit_func_p6(xarr, fit_result["a6"].median_p[0], fit_result["b6"].median[0], fit_result["c6"].median[0], fit_result["d6"].median[0])
#     yarr_6m = fit_func_p6(xarr, fit_result["a6"].median_m[0], fit_result["b6"].median[0], fit_result["c6"].median[0], fit_result["d6"].median[0])
#     print("a_0(P²): %f + %f - %f"%(fit_result["a_0_2"].median[0],fit_result["a_0_2"].ep[0],fit_result["a_0_2"].em[0]))
#     print("a_0(P⁴): %f + %f - %f"%(fit_result["a_0_4"].median[0],fit_result["a_0_4"].ep[0],fit_result["a_0_4"].em[0]))
#     print("a_0(P⁶): %f + %f - %f"%(fit_result["a_0_6"].median[0],fit_result["a_0_6"].ep[0],fit_result["a_0_6"].em[0]))
#     plt.plot(xarr, yarr_2, color = color_arr[0], label = "$\\mathcal{O}(P^2)$")#: a_0 = %1.3f +  %1.3f -  %1.3f$"%(fit_result["a_0_2"].median[0],fit_result["a_0_2"].ep[0],fit_result["a_0_2"].em[0]))
#     plt.fill_between(x=xarr, y1=yarr_2m, y2=yarr_2p, color = color_arr[0], alpha = 0.3)
#     plt.plot(xarr, yarr_4, color = color_arr[1], label = "$\\mathcal{O}(P^4)$")#: a_0 = %1.3f +  %1.3f -  %1.3f$"%(fit_result["a_0_4"].median[0],fit_result["a_0_4"].ep[0],fit_result["a_0_4"].em[0]))
#     plt.fill_between(x=xarr, y1=yarr_4m, y2=yarr_4p, color = color_arr[1], alpha = 0.3)
#     # plt.plot(xarr, yarr_6, color = color_arr[2], label = "$\\mathcal{O}(P^6): a_0 = %1.3f +  %1.3f -  %1.3f$"%(fit_result["a_0_6"].median[0],fit_result["a_0_6"].ep[0],fit_result["a_0_6"].em[0]))
#     # plt.fill_between(x=xarr, y1=yarr_6p, y2=yarr_6m, color = color_arr[2], alpha = 0.3)
#     # plt.plot(xarr, yarr_2_fixed, color = color_arr[3], label = "$\\mathcal{O}(P^2)$ (fixed)")
#     # plt.fill_between(x=xarr, y1=yarr_2m_fixed, y2=yarr_2p_fixed, color = color_arr[3], alpha = 0.3)
#     plt.plot(xarr, yarr_4_fixed, color = color_arr[4], label = "$\\mathcal{O}(P^4)$ (fixed)")
#     plt.fill_between(x=xarr, y1=yarr_4m_fixed, y2=yarr_4p_fixed, color = color_arr[4], alpha = 0.3)

# def fit_func_tan2(P2, a, b):
#     return 1/(a+b*P2)
# def fit_func_tan2_fixed_400(P2, b):
#     return 1/(-1/400+b*P2)
# def fit_func_tan4(P2, a, b, c):
#     return 1/(a+b*P2+c*P2*P2)
# def fit_func_tan4_fixed_400(P2, b, c):
#     return 1/(-1/400+b*P2+c*P2*P2)
# def fit_func_tan6(P2, a, b, c, d):
#     return 1/(a+b*P2+c*P2*P2+d*P2*P2*P2)

# def plot_fits_tan_PS(P_2, fit_result):
#     xarr = np.logspace(np.log(min(P_2)[0]*0.9),np.log(max(P_2)[0]*1.1), 500)
#     yarr_2 = fit_func_tan2(xarr, fit_result["a2"].median[0], fit_result["b2"].median[0])
#     yarr_2p = fit_func_tan2(xarr, fit_result["a2"].median_p[0], fit_result["b2"].median[0])
#     yarr_2m = fit_func_tan2(xarr, fit_result["a2"].median_m[0], fit_result["b2"].median[0])
#     yarr_2_fixed = fit_func_tan2_fixed_400(xarr, fit_result["b2_fixed"].median[0])
#     yarr_2p_fixed = fit_func_tan2_fixed_400(xarr, fit_result["b2_fixed"].median_p[0])
#     yarr_2m_fixed = fit_func_tan2_fixed_400(xarr, fit_result["b2_fixed"].median_m[0])
#     yarr_4 = fit_func_tan4(xarr, fit_result["a4"].median[0], fit_result["b4"].median[0], fit_result["c4"].median[0])
#     yarr_4p = fit_func_tan4(xarr, fit_result["a4"].median_p[0], fit_result["b4"].median[0], fit_result["c4"].median[0])
#     yarr_4m = fit_func_tan4(xarr, fit_result["a4"].median_m[0], fit_result["b4"].median[0], fit_result["c4"].median[0])
#     yarr_4_fixed = fit_func_tan4_fixed_400(xarr, fit_result["b4_fixed"].median[0], fit_result["c4_fixed"].median[0])
#     yarr_4p_fixed = fit_func_tan4_fixed_400(xarr, fit_result["b4_fixed"].median_p[0], fit_result["c4_fixed"].median[0])
#     yarr_4m_fixed = fit_func_tan4_fixed_400(xarr, fit_result["b4_fixed"].median_m[0], fit_result["c4_fixed"].median[0])
#     yarr_6 = fit_func_tan6(xarr, fit_result["a6"].median[0], fit_result["b6"].median[0], fit_result["c6"].median[0], fit_result["d6"].median[0])
#     yarr_6p = fit_func_tan6(xarr, fit_result["a6"].median_p[0], fit_result["b6"].median[0], fit_result["c6"].median[0], fit_result["d6"].median[0])
#     yarr_6m = fit_func_tan6(xarr, fit_result["a6"].median_m[0], fit_result["b6"].median[0], fit_result["c6"].median[0], fit_result["d6"].median[0])
#     print("a_0(P²): %f + %f - %f"%(fit_result["a_0_2"].median[0],fit_result["a_0_2"].ep[0],fit_result["a_0_2"].em[0]))
#     print("a_0(P⁴): %f + %f - %f"%(fit_result["a_0_4"].median[0],fit_result["a_0_4"].ep[0],fit_result["a_0_4"].em[0]))
#     print("a_0(P⁶): %f + %f - %f"%(fit_result["a_0_6"].median[0],fit_result["a_0_6"].ep[0],fit_result["a_0_6"].em[0]))
#     plt.plot(xarr, yarr_2, color = color_arr[0], label = "$\\mathcal{O}(P^2)$") #: a_0 = %1.3f +  %1.3f -  %1.3f$"%(fit_result["a_0_2"].median[0],fit_result["a_0_2"].ep[0],fit_result["a_0_2"].em[0]))
#     plt.fill_between(x=xarr, y1=yarr_2m, y2=yarr_2p, color = color_arr[0], alpha = 0.3)
#     # plt.plot(xarr, yarr_2_fixed, color = color_arr[3], label = "$\\mathcal{O}(P^2)$ (fixed)")
#     # plt.fill_between(x=xarr, y1=yarr_2m_fixed, y2=yarr_2p_fixed, color = color_arr[3], alpha = 0.3)
#     plt.plot(xarr, yarr_4, color = color_arr[1], label = "$\\mathcal{O}(P^4)$") #: a_0 = %1.3f +  %1.3f -  %1.3f$"%(fit_result["a_0_4"].median[0],fit_result["a_0_4"].ep[0],fit_result["a_0_4"].em[0]))
#     plt.fill_between(x=xarr, y1=yarr_4m, y2=yarr_4p, color = color_arr[1], alpha = 0.3)
#     plt.plot(xarr, yarr_4_fixed, color = color_arr[4], label = "$\\mathcal{O}(P^4)$ (fixed)")
#     plt.fill_between(x=xarr, y1=yarr_4m_fixed, y2=yarr_4p_fixed, color = color_arr[4], alpha = 0.3)
#     # plt.plot(xarr, yarr_6, color = color_arr[2], label = "$\\mathcal{O}(P^6): a_0 = %1.3f +  %1.3f -  %1.3f$"%(fit_result["a_0_6"].median[0],fit_result["a_0_6"].ep[0],fit_result["a_0_6"].em[0]))
#     # plt.fill_between(x=xarr, y1=yarr_6p, y2=yarr_6m, color = color_arr[2], alpha = 0.3)

# def plot_P_cot_PS(fit = False, zoom = [0.008,10,-0.11,0.2], title = "", save_app = "", log = ""):
#     results = get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot")
#     # results = get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot_SU3")
#     # fit_result = get_result_file("phase_shift_fit_P_cot_PS")

#     P_cot_PS, P_cot_PS_err, P_2, P_2_err, rho_pi = [[],[],[],[],[]]
#     for i in range(len(results)):
#         rho_pi.append(results[i]["infos/m_rho_pi"])
#         P_2.append(results[i]["P_2"].median)
#         print("E_pipi: ", results[i]["E_pipi"].median)
#         print("tan_PS: ", results[i]["tan_PS"].median)
#         print("cot_PS: ", results[i]["cot_PS"].median)
#         print("P: ", results[i]["P"].median)
#         print("P_2: ", results[i]["P_2"].median)
#         print("m_inf: ", results[i]["mass_Goldstone"].median)
#         print("N_L: ", results[i]["N_L"].median)
#         print("\n")
#         P_2_err.append([results[i]["P_2"].ep,results[i]["P_2"].em])
#         P_cot_PS.append(results[i]["cot_PS"].median[0])
#         P_cot_PS_err.append([results[i]["cot_PS"].ep,results[i]["cot_PS"].em])
    
#     norm = matplotlib.colors.Normalize(vmin=min(rho_pi), vmax=max(rho_pi), clip=True)
#     mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
#     time_color = np.array([(mapper.to_rgba(v)) for v in rho_pi])
#     plt.colorbar(mappable=mapper, label = "$\\frac{m_{\\rho}}{m_{\\pi}}$")
#     for i in range(len(results)):
#         # plt.errorbar(x=P_2[i], xerr = P_2_err[i], y=P_cot_PS[i], yerr=P_cot_PS_err[i], color = time_color[i], ls = "", capsize=5, markersize=10)
#         plt.scatter(x=P_2[i], y=P_cot_PS[i], color = time_color[i])
#     if fit:
#         plot_fits_P_cot_PS(P_2, fit_result)


#     # plt.axvline(0.4, color = "black", ls="--")
#     log_txt = ""
#     if log:
#         plt.xscale("log")
#         log_txt = "_log"
#     plt.ylabel("$\\frac{Pcot(\delta)}{m_{\pi}}$")
#     plt.xlabel("$\\frac{P^2}{m_{\pi}^2}$")
#     plt.legend(fontsize = "x-small")
#     plt.title(title)
#     plt.grid()
#     zoom_txt = ""
#     if zoom == plot_tan_PS.__defaults__[0]:
#         zoom_txt = "_zoom"
#     if zoom != False:
#         plt.xlim((zoom[0], zoom[1]))
#         plt.ylim((zoom[2], zoom[3]))
#     plt.savefig("plots/P_cot_PS"+log_txt+save_app+zoom_txt+".pdf")
#     plt.show()
#     plt.clf()

# def plot_tan_PS(fit = False, zoom = [0.008,10,-400,100], title = "", save_app = "", log = False):
#     results = get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot")
#     # results = get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot_CDR")
#     # results = get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot_SU3")
#     # fit_result = get_result_file("phase_shift_fit_P_cot_PS")

#     tan_PS, tan_PS_err, P_2, P_2_err, rho_pi = [[],[],[],[],[]]
#     for i in range(len(results)):
#         rho_pi.append(results[i]["infos/m_rho_pi"])
#         P_2.append(results[i]["P"].median[0])
#         P_2_err.append([results[i]["P"].ep,results[i]["P"].em])
#         tan_PS.append(results[i]["tan_PS"].median[0])
#         tan_PS_err.append([results[i]["tan_PS"].ep,results[i]["tan_PS"].em])
#     norm = matplotlib.colors.Normalize(vmin=min(rho_pi), vmax=max(rho_pi), clip=True)
#     mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
#     time_color = np.array([(mapper.to_rgba(v)) for v in rho_pi])
#     plt.colorbar(mappable=mapper, label = "$\\frac{m_{\\rho}}{m_{\\pi}}$")
#     for i in range(len(results)):
#         print(P_2[i], tan_PS[i])
#         # plt.errorbar(x=P_2[i], xerr = P_2_err[i], y=tan_PS[i], yerr=tan_PS_err[i], color = time_color[i], ls = "", capsize=5, markersize=10)
#         plt.scatter(x=P_2[i], y=tan_PS[i], color = time_color[i])

#     if fit:
#         plot_fits_tan_PS(P_2, fit_result)

#     plt.title(title)
#     log_txt = ""
#     if log:
#         plt.xscale("log")
#         log_txt = "_log"
#     plt.ylabel("$tan(\delta)\\frac{m_{\pi}}{P}$")
#     plt.xlabel("$\\frac{P^2}{m_{\pi}^2}$")
#     plt.legend(fontsize = "x-small")
#     plt.grid()
#     zoom_txt = ""
#     if zoom == plot_tan_PS.__defaults__[0]:
#         zoom_txt = "_zoom"
#     plt.savefig("plots/tan_PS"+log_txt+save_app+zoom_txt+".pdf")
#     plt.show()
#     plt.clf()






# def main():
#     # create_all_filenames()
#     # plot_P_cot_PS(fit = False, zoom = False)
#     # plot_tan_PS(fit = False, zoom = False, log = True)
#     # plot_P_cot_PS(fit = False, save_app = "_SU3", zoom = False)
#     # plot_tan_PS(fit = False, zoom=[0,0.5,-4,1], save_app = "_SU3_tan_P", log=True)
#     # plot_tan_PS(fit = False, save_app = "_SU3", log=True)

#     plot_tan_PS(fit = False, zoom = False, save_app="_SP4", log = False)
#     # plot_P_cot_PS(fit = False, zoom = False, save_app="_SP4", log = False)
#     # plot_tan_PS(fit = False, zoom=[0,0.5,-4,1], save_app = "_SU3_tan", log=True)


#     # plot_P_cot_PS(fit = True, title="P cot vs P², Data don't show y-axis crossing", zoom=[-0.04,0.2,-0.06,0.01],save_app="1")
#     # plot_P_cot_PS(fit = True, title="P cot vs P², Data don't show y-axis crossing", zoom=[-0.05,0.6,-0.15,0.01],save_app="2")
#     # plot_tan_PS(fit = True, title="Data is consistent with slope getting more and more negative", zoom=[-0.02,1e-1,-500,50])
#     # plot_tan_PS(fit = True, title="The fitting function demands a constant emerging for small P²\nThis is not happening in the data. Set to -400 here", zoom=[1e-4,1e1,-500,100],log="log")
#     # plot_P_cot_PS(fit = True, title="In this plot it looks like 1/a_0 is consistent with 0", save_app="fixed", zoom=[1e-4,1e1,-0.15,0.05],log="log")
    
#     # plot_P_cot_PS(fit = True, title="Der Vergleich mit der fixed Methode im cot-plot", save_app="fixed", zoom=[1e-4,1e1,-0.15,0.05])
#     # plot_P_cot_PS(fit = True, title="Hier sieht man, dass höhere Ordnungen den Fit nicht verbessern", zoom=[1e-4,1e1,-0.15,0.05],log="log")
#     # plot_zeta()
#     # plot_q_P()
#     # plot_dispersion_relation()


# def print_phase_shift_data():
#     with open("output/phase_shift_data.dat", "w") as file:
#         with open("output/phase_shift_data_mathematica.dat", "w") as file_mathematica:
#             file.write("P_2\tP_2_err\P_cot\tP_cot_err\ttan/P\ttan/P_err\tm_rho/m_pi\n")
#             file_mathematica.write("i\tP_2\tP_2_err\P_cot\tP_cot_err\ttan/P\ttan/P_err\tm_rho/m_pi\n{")
#             i = 1
#             for result in get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot"):
#                 rho_pi = result["infos/m_rho_pi"]
#                 file.write("%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\n"%(result["P_2_prime"].median[0],result["P_2_prime"].e[0],result["P_cot_PS_prime"].median[0],result["P_cot_PS_prime"].e[0],result["tan_PS/P_prime"].median[0],result["P_cot_PS_prime"].e[0],rho_pi))

#                 rho_pi = result["infos/m_rho_pi"]
#                 file_mathematica.write("{%i,%.30f,%.30f,%.30f,%.30f,%.30f,%.30f,%.30f},"%(i,result["P_2_prime"].median[0],result["P_2_prime"].e[0],result["P_cot_PS_prime"].median[0],result["P_cot_PS_prime"].e[0],result["tan_PS/P_prime"].median[0],result["P_cot_PS_prime"].e[0],rho_pi))                                  # DONT forget to erase last comma
#                 i = i+1
#             file_mathematica.write("}")

# def print_energy_level_data():
#     # with open("output/energy_level_data.dat", "w") as file:
#     with open("output/energy_level_data_short.dat", "w") as file:
#         # file.write("gauge_group\tbeta\tm_1\tm_2\tN_T\tN_L\tlevel\tm_rho_pi\tm_pi_inf\tm_pi_inf_err\tE\t\tE_err\n")
#         for result in get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_all"):
#             for i in range(len(result["E"].median)):
#                 file.write("%s\t\t"%result["infos/gauge_group"])
#                 file.write("%1.3f\t"%result["infos/beta"])
#                 file.write("%1.3f\t"%result["infos/m_1"])
#                 file.write("%1.3f\t"%result["infos/m_2"])
#                 file.write("%i\t"%result["N_Ts"].median[i])
#                 file.write("%i\t"%result["N_Ls"].median[i])
#                 file.write("%i\t"%result["infos/level"])
#                 file.write("%e\t"%result["infos/m_rho_pi"])
#                 file.write("%e\t"%result["infos/m_pi_inf"])
#                 file.write("%e\t"%result["infos/m_pi_inf_err"])
#                 file.write("%e\t"%result["E"].median[i])
#                 file.write("%e\t"%result["E"].e[i])
#                 file.write("\n")
#             # file.write("\n")
            
#                 # file.write("\n___________________________________________________________________________________________________\n")
                
#                 # print(result.keys(),"\n")
#                 # file.write("P_2\tP_2_err\P_cot\tP_cot_err\ttan/P\ttan/P_err\tm_rho/m_pi\n")

    # print_phase_shift_data()
    # print_energy_level_data()












    # for resultfile in resultfile_list:
    #     phase_shift = errcl.measurement(resultfile)
    #     phase_shift.read_from_HDF()
    #     info = phase_shift.infos
    #     pi_energy_lev = errcl.measurement("energy_levels_Fabian_Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i_pi"%(info["isospin_channel"],info["gauge_group"],info["beta"],info["m_1"],info["m_2"],info["N_T"],info["N_L"]))
    #     pi_energy_lev.read_from_HDF()
    #     rho_energy_lev = errcl.measurement("energy_levels_Fabian_Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i_rho"%(info["isospin_channel"],info["gauge_group"],info["beta"],info["m_1"],info["m_2"],info["N_T"],info["N_L"]))
    #     rho_energy_lev.read_from_HDF()
    #     # pi_energy_lev.print_everything()
    #     rho_pi.append(rho_energy_lev.results["E"].median[0]/pi_energy_lev.results["E"].median[0])
    #     infos = phase_shift.infos
    #     N_T = infos["N_T"]
    #     N_L = infos["N_L"]
    #     info_str = "Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(infos["isospin_channel"],infos["gauge_group"],infos["beta"],infos["m_1"],infos["m_2"])
    #     ensemble_str = "beta=%1.3f, m1/2=%1.3f"%(infos["beta"],infos["m_1"])
    #     y.append(phase_shift.results["P_cot_PS_prime"].median)
    #     yerr.append([phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em])
    #     x.append(phase_shift.results["P_2_prime"].median)
    #     xerr.append([phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em])

    # phase_shift_fit = errcl.measurement("phase_shift_fit_P_cot_PS")
    # phase_shift_fit.read_from_HDF()
    # xarr = np.logspace(np.log(min(x)[0]*0.9),np.log(max(x)[0]*1.1), 500)
    # # xarr = np.linspace(min(min(x))*0.9, max(max(x))*1.1, 500)
    # yarr_2 = fit_func_p2(xarr, phase_shift_fit.results["a2"].median[0], phase_shift_fit.results["b2"].median[0])
    # yarr_2p = fit_func_p2(xarr, phase_shift_fit.results["a2"].median_p[0], phase_shift_fit.results["b2"].median[0])
    # yarr_2m = fit_func_p2(xarr, phase_shift_fit.results["a2"].median_m[0], phase_shift_fit.results["b2"].median[0])
    # yarr_4 = fit_func_p4(xarr, phase_shift_fit.results["a4"].median[0], phase_shift_fit.results["b4"].median[0], phase_shift_fit.results["c4"].median[0])
    # yarr_4p = fit_func_p4(xarr, phase_shift_fit.results["a4"].median_p[0], phase_shift_fit.results["b4"].median[0], phase_shift_fit.results["c4"].median[0])
    # yarr_4m = fit_func_p4(xarr, phase_shift_fit.results["a4"].median_m[0], phase_shift_fit.results["b4"].median[0], phase_shift_fit.results["c4"].median[0])
    # yarr_6 = fit_func_p6(xarr, phase_shift_fit.results["a6"].median[0], phase_shift_fit.results["b6"].median[0], phase_shift_fit.results["c6"].median[0], phase_shift_fit.results["d6"].median[0])
    # yarr_6p = fit_func_p6(xarr, phase_shift_fit.results["a6"].median_p[0], phase_shift_fit.results["b6"].median[0], phase_shift_fit.results["c6"].median[0], phase_shift_fit.results["d6"].median[0])
    # yarr_6m = fit_func_p6(xarr, phase_shift_fit.results["a6"].median_m[0], phase_shift_fit.results["b6"].median[0], phase_shift_fit.results["c6"].median[0], phase_shift_fit.results["d6"].median[0])

    # print("a_0(P²): %f + %f - %f"%(phase_shift_fit.results["a_0_2"].median[0],phase_shift_fit.results["a_0_2"].ep[0],phase_shift_fit.results["a_0_2"].em[0]))
    # print("a_0(P⁴): %f + %f - %f"%(phase_shift_fit.results["a_0_4"].median[0],phase_shift_fit.results["a_0_4"].ep[0],phase_shift_fit.results["a_0_4"].em[0]))
    # print("a_0(P⁶): %f + %f - %f"%(phase_shift_fit.results["a_0_6"].median[0],phase_shift_fit.results["a_0_6"].ep[0],phase_shift_fit.results["a_0_6"].em[0]))

    # fig = plt.figure()
    # ax = fig.add_subplot(111)

    # plt.plot(xarr, yarr_2, color = color_arr[0], label = "$\\mathcal{O}(P^2): a_0 = %1.3f +  %1.3f -  %1.3f$"%(phase_shift_fit.results["a_0_2"].median[0],phase_shift_fit.results["a_0_2"].ep[0],phase_shift_fit.results["a_0_2"].em[0]))
    # plt.fill_between(x=xarr, y1=yarr_2m, y2=yarr_2p, color = color_arr[0], alpha = 0.3)
    # plt.plot(xarr, yarr_4, color = color_arr[1], label = "$\\mathcal{O}(P^4): a_0 = %1.3f +  %1.3f -  %1.3f$"%(phase_shift_fit.results["a_0_4"].median[0],phase_shift_fit.results["a_0_4"].ep[0],phase_shift_fit.results["a_0_4"].em[0]))
    # plt.fill_between(x=xarr, y1=yarr_4m, y2=yarr_4p, color = color_arr[1], alpha = 0.3)
    # plt.plot(xarr, yarr_6, color = color_arr[2], label = "$\\mathcal{O}(P^6): a_0 = %1.3f +  %1.3f -  %1.3f$"%(phase_shift_fit.results["a_0_6"].median[0],phase_shift_fit.results["a_0_6"].ep[0],phase_shift_fit.results["a_0_6"].em[0]))
    # plt.fill_between(x=xarr, y1=yarr_6p, y2=yarr_6m, color = color_arr[2], alpha = 0.3)

    # norm = matplotlib.colors.Normalize(vmin=min(rho_pi), vmax=max(rho_pi), clip=True)
    # mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    # time_color = np.array([(mapper.to_rgba(v)) for v in rho_pi])
    # plt.colorbar(mappable=mapper, label = "$\\frac{m_{\\rho}}{m_{\\pi}}$")
    # for i in range(len(yerr)):
    #     plt.errorbar(x=x[i], xerr = xerr[i], y=y[i], yerr=yerr[i], color = time_color[i], ls = "", capsize=5, markersize=10)
    # plt.xlim((zoom[0], zoom[1]))
    # plt.ylim((zoom[2], zoom[3]))

    # plt.axvline(0.4, color = "black", ls="--")
    # plt.xscale("log")
    # plt.ylabel("$\\frac{Pcot(\delta)}{m_{\pi}}$")
    # plt.xlabel("$\\frac{P^2}{m_{\pi}^2}$")
    # plt.legend(fontsize = "x-small")
    # plt.grid()
    # # plt.title(info[7])
    # if zoom == plot_P_cot_PS_fit.__defaults__[0]:
    #     plt.savefig("plots/P_cot_PS_fit.pdf")
    # else:
    #     plt.savefig("plots/P_cot_PS_fit_zoom.pdf")
    # plt.show()
    # plt.clf()





















































# def plot_phase_shift():
#     ensemble_list = []
#     color_ind = -1


#     y = []
#     yerr = []
#     x = []
#     xerr = []
#     rho_pi = []

#     # resultfile_list = get_result_files("phase_shift")
#     # resultfile_list = get_result_files("phase_shift_const")
#     resultfile_list = get_result_files("phase_shift_Fabian")
#     for resultfile in resultfile_list:
#         phase_shift = errcl.measurement(resultfile)
#         phase_shift.read_from_HDF()
#         info = phase_shift.infos
#         pi_energy_lev = errcl.measurement("energy_levels_Fabian_Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i_pi"%(info["isospin_channel"],info["gauge_group"],info["beta"],info["m_1"],info["m_2"],info["N_T"],info["N_L"]))
#         pi_energy_lev.read_from_HDF()
#         rho_energy_lev = errcl.measurement("energy_levels_Fabian_Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i_rho"%(info["isospin_channel"],info["gauge_group"],info["beta"],info["m_1"],info["m_2"],info["N_T"],info["N_L"]))
#         rho_energy_lev.read_from_HDF()
#         pi_energy_lev.print_everything()
#         rho_pi.append(rho_energy_lev.results["E"].median[0]/pi_energy_lev.results["E"].median[0])
#         infos = phase_shift.infos
#         N_T = infos["N_T"]
#         N_L = infos["N_L"]
#         info_str = "Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(infos["isospin_channel"],infos["gauge_group"],infos["beta"],infos["m_1"],infos["m_2"])
#         ensemble_str = "beta=%1.3f, m1/2=%1.3f"%(infos["beta"],infos["m_1"])
#         y.append(phase_shift.results["P_cot_PS_prime"].median)
#         yerr.append([phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em])
#         x.append(phase_shift.results["P_2_prime"].median)
#         xerr.append([phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em])

#     norm = matplotlib.colors.Normalize(vmin=min(rho_pi), vmax=max(rho_pi), clip=True)
#     mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
#     time_color = np.array([(mapper.to_rgba(v)) for v in rho_pi])
#     for i in range(len(yerr)):
#         plt.errorbar(x=x[i], xerr = xerr[i], y=y[i], yerr=yerr[i], color = time_color[i], ls = "", capsize=5, markersize=10)
        
#     plt.xscale("log")
#     plt.ylabel("$\\frac{Pcot(\delta)}{m_{\pi}}$")
#     plt.xlabel("$\\frac{P^2}{m_{\pi}}$")
#     plt.legend(fontsize = "xx-small")
#     plt.xlim((0.008, 10))
#     plt.ylim((-0.11,0.2))
#     plt.grid()
#     # plt.title(info[7])
#     plt.savefig("plots/P_cot_PS_.pdf")
#     plt.show()
#     plt.clf()

# def fit_func_p2(P2, a, b):
#     return a + b*P2
# def fit_func_p4(P2, a, b, c):
#     return a + b*P2 + c*P2*P2
# def fit_func_p6(P2, a, b, c, d):
#     return a + b*P2 + c*P2*P2 + d*P2*P2*P2

# def plot_P_cot_PS_fit(zoom = [0.008,10,-0.11,0.2]):
#     ensemble_list = []
#     color_ind = -1

#     y = []
#     yerr = []
#     x = []
#     xerr = []
#     rho_pi = []

#     # resultfile_list = get_result_files("phase_shift")
#     # resultfile_list = get_result_files("phase_shift_const")
#     resultfile_list = get_result_files("phase_shift_Fabian")
#     for resultfile in resultfile_list:
#         phase_shift = errcl.measurement(resultfile)
#         phase_shift.read_from_HDF()
#         info = phase_shift.infos
#         pi_energy_lev = errcl.measurement("energy_levels_Fabian_Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i_pi"%(info["isospin_channel"],info["gauge_group"],info["beta"],info["m_1"],info["m_2"],info["N_T"],info["N_L"]))
#         pi_energy_lev.read_from_HDF()
#         rho_energy_lev = errcl.measurement("energy_levels_Fabian_Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i_rho"%(info["isospin_channel"],info["gauge_group"],info["beta"],info["m_1"],info["m_2"],info["N_T"],info["N_L"]))
#         rho_energy_lev.read_from_HDF()
#         # pi_energy_lev.print_everything()
#         rho_pi.append(rho_energy_lev.results["E"].median[0]/pi_energy_lev.results["E"].median[0])
#         infos = phase_shift.infos
#         N_T = infos["N_T"]
#         N_L = infos["N_L"]
#         info_str = "Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(infos["isospin_channel"],infos["gauge_group"],infos["beta"],infos["m_1"],infos["m_2"])
#         ensemble_str = "beta=%1.3f, m1/2=%1.3f"%(infos["beta"],infos["m_1"])
#         y.append(phase_shift.results["P_cot_PS_prime"].median)
#         yerr.append([phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em])
#         x.append(phase_shift.results["P_2_prime"].median)
#         xerr.append([phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em])

#     phase_shift_fit = errcl.measurement("phase_shift_fit_P_cot_PS")
#     phase_shift_fit.read_from_HDF()
#     xarr = np.logspace(np.log(min(x)[0]*0.9),np.log(max(x)[0]*1.1), 500)
#     # xarr = np.linspace(min(min(x))*0.9, max(max(x))*1.1, 500)
#     yarr_2 = fit_func_p2(xarr, phase_shift_fit.results["a2"].median[0], phase_shift_fit.results["b2"].median[0])
#     yarr_2p = fit_func_p2(xarr, phase_shift_fit.results["a2"].median_p[0], phase_shift_fit.results["b2"].median[0])
#     yarr_2m = fit_func_p2(xarr, phase_shift_fit.results["a2"].median_m[0], phase_shift_fit.results["b2"].median[0])
#     yarr_4 = fit_func_p4(xarr, phase_shift_fit.results["a4"].median[0], phase_shift_fit.results["b4"].median[0], phase_shift_fit.results["c4"].median[0])
#     yarr_4p = fit_func_p4(xarr, phase_shift_fit.results["a4"].median_p[0], phase_shift_fit.results["b4"].median[0], phase_shift_fit.results["c4"].median[0])
#     yarr_4m = fit_func_p4(xarr, phase_shift_fit.results["a4"].median_m[0], phase_shift_fit.results["b4"].median[0], phase_shift_fit.results["c4"].median[0])
#     yarr_6 = fit_func_p6(xarr, phase_shift_fit.results["a6"].median[0], phase_shift_fit.results["b6"].median[0], phase_shift_fit.results["c6"].median[0], phase_shift_fit.results["d6"].median[0])
#     yarr_6p = fit_func_p6(xarr, phase_shift_fit.results["a6"].median_p[0], phase_shift_fit.results["b6"].median[0], phase_shift_fit.results["c6"].median[0], phase_shift_fit.results["d6"].median[0])
#     yarr_6m = fit_func_p6(xarr, phase_shift_fit.results["a6"].median_m[0], phase_shift_fit.results["b6"].median[0], phase_shift_fit.results["c6"].median[0], phase_shift_fit.results["d6"].median[0])

#     print("a_0(P²): %f + %f - %f"%(phase_shift_fit.results["a_0_2"].median[0],phase_shift_fit.results["a_0_2"].ep[0],phase_shift_fit.results["a_0_2"].em[0]))
#     print("a_0(P⁴): %f + %f - %f"%(phase_shift_fit.results["a_0_4"].median[0],phase_shift_fit.results["a_0_4"].ep[0],phase_shift_fit.results["a_0_4"].em[0]))
#     print("a_0(P⁶): %f + %f - %f"%(phase_shift_fit.results["a_0_6"].median[0],phase_shift_fit.results["a_0_6"].ep[0],phase_shift_fit.results["a_0_6"].em[0]))

#     fig = plt.figure()
#     ax = fig.add_subplot(111)

#     plt.plot(xarr, yarr_2, color = color_arr[0], label = "$\\mathcal{O}(P^2): a_0 = %1.3f +  %1.3f -  %1.3f$"%(phase_shift_fit.results["a_0_2"].median[0],phase_shift_fit.results["a_0_2"].ep[0],phase_shift_fit.results["a_0_2"].em[0]))
#     plt.fill_between(x=xarr, y1=yarr_2m, y2=yarr_2p, color = color_arr[0], alpha = 0.3)
#     plt.plot(xarr, yarr_4, color = color_arr[1], label = "$\\mathcal{O}(P^4): a_0 = %1.3f +  %1.3f -  %1.3f$"%(phase_shift_fit.results["a_0_4"].median[0],phase_shift_fit.results["a_0_4"].ep[0],phase_shift_fit.results["a_0_4"].em[0]))
#     plt.fill_between(x=xarr, y1=yarr_4m, y2=yarr_4p, color = color_arr[1], alpha = 0.3)
#     plt.plot(xarr, yarr_6, color = color_arr[2], label = "$\\mathcal{O}(P^6): a_0 = %1.3f +  %1.3f -  %1.3f$"%(phase_shift_fit.results["a_0_6"].median[0],phase_shift_fit.results["a_0_6"].ep[0],phase_shift_fit.results["a_0_6"].em[0]))
#     plt.fill_between(x=xarr, y1=yarr_6p, y2=yarr_6m, color = color_arr[2], alpha = 0.3)

#     norm = matplotlib.colors.Normalize(vmin=min(rho_pi), vmax=max(rho_pi), clip=True)
#     mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
#     time_color = np.array([(mapper.to_rgba(v)) for v in rho_pi])
#     plt.colorbar(mappable=mapper, label = "$\\frac{m_{\\rho}}{m_{\\pi}}$")
#     for i in range(len(yerr)):
#         plt.errorbar(x=x[i], xerr = xerr[i], y=y[i], yerr=yerr[i], color = time_color[i], ls = "", capsize=5, markersize=10)
#     plt.xlim((zoom[0], zoom[1]))
#     plt.ylim((zoom[2], zoom[3]))

#     plt.axvline(0.4, color = "black", ls="--")
#     plt.xscale("log")
#     plt.ylabel("$\\frac{Pcot(\delta)}{m_{\pi}}$")
#     plt.xlabel("$\\frac{P^2}{m_{\pi}^2}$")
#     plt.legend(fontsize = "x-small")
#     plt.grid()
#     # plt.title(info[7])
#     if zoom == plot_P_cot_PS_fit.__defaults__[0]:
#         plt.savefig("plots/P_cot_PS_fit.pdf")
#     else:
#         plt.savefig("plots/P_cot_PS_fit_zoom.pdf")
#     plt.show()
#     plt.clf()





# def main():
#     # plot_basic_analysis()
#     # plot_energy_levels()
#     # plot_energy_levels_Fabian()
#     # plot_infinite_volume()
#     # compare_energy_levels()
#     # plot_phase_shift()
#     plot_P_cot_PS_fit()
#     plot_P_cot_PS_fit(zoom=[0.01,0.6,-0.11, 0.01])
#     print_phase_shift_data()


# if __name__ == "__main__":
#     main()