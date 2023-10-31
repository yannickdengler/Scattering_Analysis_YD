import os
import matplotlib.pyplot as plt
import read_HDF5_logfile as HDF_log
import error_classes as errcl 
import numpy as np
import matplotlib
from matplotlib import cm

color_arr = ["blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick"]

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
        results.append(get_result_file(filelist[i]))
    return results

def get_result_file(file):
    dic = {}
    meas = errcl.measurement(file)
    meas.read_from_HDF()
    for key in meas.result_names:
        # print(key)
        dic[key] = meas.results[key]
    for key in meas.infos.keys():
        dic["infos/"+key] = meas.infos[key]
    return dic

def plot_zeta():
    import generalizedzeta as gz
    results = get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot")
    
    zeta, zeta_err, q_2, q_2_err, rho_pi = [[],[],[],[],[]]
    for i in range(len(results)):
        rho_pi.append(results[i]["infos/m_rho_pi"])
        zeta.append(results[i]["Zeta"].median)
        zeta_err.append([results[i]["Zeta"].ep,results[i]["Zeta"].em])
        q_2.append(results[i]["q_2"].median[0])
        q_2_err.append([results[i]["q_2"].ep,results[i]["q_2"].em])
    # xarr = np.linspace(min(q_2)*0.9, max(q_2)*1.1, 500)
    xarr = np.linspace(-1, 1.5, 500)
    yarr = []
    for x in xarr:
        yarr.append(gz.Zeta(x))
    plt.plot(xarr, yarr)

    norm = matplotlib.colors.Normalize(vmin=min(rho_pi), vmax=max(rho_pi), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    time_color = np.array([(mapper.to_rgba(v)) for v in rho_pi])
    plt.colorbar(mappable=mapper, label = "$\\frac{m_{\\rho}}{m_{\\pi}}$")
    for i in range(len(results)):
        plt.errorbar(x=q_2[i], xerr = q_2_err[i], y=zeta[i], yerr=zeta_err[i], color = time_color[i], ls = "", capsize=5, markersize=10)
    plt.grid()
    plt.xlim((-0.5,1.5))
    plt.ylim((-55,10))
    plt.savefig("Zeta.pdf")
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



def fit_func_p2(P2, a, b):
    return a + b*P2
def fit_func_p4(P2, a, b, c):
    return a + b*P2 + c*P2*P2
def fit_func_p6(P2, a, b, c, d):
    return a + b*P2 + c*P2*P2 + d*P2*P2*P2
def fit_func_p2_fixed_400(P2, b):
    return -1/400 + b*P2
def fit_func_p4_fixed_400(P2, b, c):
    return -1/400 + b*P2 + c*P2*P2

def plot_fits_P_cot_PS(P_2, fit_result):
    xarr = np.logspace(np.log(min(P_2)[0]*0.9),np.log(max(P_2)[0]*1.1), 500)
    yarr_2 = fit_func_p2(xarr, fit_result["a2"].median[0], fit_result["b2"].median[0])
    yarr_2p = fit_func_p2(xarr, fit_result["a2"].median_p[0], fit_result["b2"].median[0])
    yarr_2m = fit_func_p2(xarr, fit_result["a2"].median_m[0], fit_result["b2"].median[0])
    yarr_2_fixed = fit_func_p2_fixed_400(xarr, fit_result["b2_fixed"].median[0])
    yarr_2p_fixed = fit_func_p2_fixed_400(xarr, fit_result["b2_fixed"].median[0])
    yarr_2m_fixed = fit_func_p2_fixed_400(xarr, fit_result["b2_fixed"].median[0])
    yarr_4 = fit_func_p4(xarr, fit_result["a4"].median[0], fit_result["b4"].median[0], fit_result["c4"].median[0])
    yarr_4p = fit_func_p4(xarr, fit_result["a4"].median_p[0], fit_result["b4"].median[0], fit_result["c4"].median[0])
    yarr_4m = fit_func_p4(xarr, fit_result["a4"].median_m[0], fit_result["b4"].median[0], fit_result["c4"].median[0])
    yarr_4_fixed = fit_func_p4_fixed_400(xarr, fit_result["b4_fixed"].median[0], fit_result["c4_fixed"].median[0])
    yarr_4p_fixed = fit_func_p4_fixed_400(xarr, fit_result["b4_fixed"].median_p[0], fit_result["c4_fixed"].median[0])
    yarr_4m_fixed = fit_func_p4_fixed_400(xarr, fit_result["b4_fixed"].median_m[0], fit_result["c4_fixed"].median[0])
    yarr_6 = fit_func_p6(xarr, fit_result["a6"].median[0], fit_result["b6"].median[0], fit_result["c6"].median[0], fit_result["d6"].median[0])
    yarr_6p = fit_func_p6(xarr, fit_result["a6"].median_p[0], fit_result["b6"].median[0], fit_result["c6"].median[0], fit_result["d6"].median[0])
    yarr_6m = fit_func_p6(xarr, fit_result["a6"].median_m[0], fit_result["b6"].median[0], fit_result["c6"].median[0], fit_result["d6"].median[0])
    print("a_0(P²): %f + %f - %f"%(fit_result["a_0_2"].median[0],fit_result["a_0_2"].ep[0],fit_result["a_0_2"].em[0]))
    print("a_0(P⁴): %f + %f - %f"%(fit_result["a_0_4"].median[0],fit_result["a_0_4"].ep[0],fit_result["a_0_4"].em[0]))
    print("a_0(P⁶): %f + %f - %f"%(fit_result["a_0_6"].median[0],fit_result["a_0_6"].ep[0],fit_result["a_0_6"].em[0]))
    plt.plot(xarr, yarr_2, color = color_arr[0], label = "$\\mathcal{O}(P^2)$")#: a_0 = %1.3f +  %1.3f -  %1.3f$"%(fit_result["a_0_2"].median[0],fit_result["a_0_2"].ep[0],fit_result["a_0_2"].em[0]))
    plt.fill_between(x=xarr, y1=yarr_2m, y2=yarr_2p, color = color_arr[0], alpha = 0.3)
    plt.plot(xarr, yarr_4, color = color_arr[1], label = "$\\mathcal{O}(P^4)$")#: a_0 = %1.3f +  %1.3f -  %1.3f$"%(fit_result["a_0_4"].median[0],fit_result["a_0_4"].ep[0],fit_result["a_0_4"].em[0]))
    plt.fill_between(x=xarr, y1=yarr_4m, y2=yarr_4p, color = color_arr[1], alpha = 0.3)
    # plt.plot(xarr, yarr_6, color = color_arr[2], label = "$\\mathcal{O}(P^6): a_0 = %1.3f +  %1.3f -  %1.3f$"%(fit_result["a_0_6"].median[0],fit_result["a_0_6"].ep[0],fit_result["a_0_6"].em[0]))
    # plt.fill_between(x=xarr, y1=yarr_6p, y2=yarr_6m, color = color_arr[2], alpha = 0.3)
    # plt.plot(xarr, yarr_2_fixed, color = color_arr[3], label = "$\\mathcal{O}(P^2)$ (fixed)")
    # plt.fill_between(x=xarr, y1=yarr_2m_fixed, y2=yarr_2p_fixed, color = color_arr[3], alpha = 0.3)
    plt.plot(xarr, yarr_4_fixed, color = color_arr[4], label = "$\\mathcal{O}(P^4)$ (fixed)")
    plt.fill_between(x=xarr, y1=yarr_4m_fixed, y2=yarr_4p_fixed, color = color_arr[4], alpha = 0.3)

def fit_func_tan2(P2, a, b):
    return 1/(a+b*P2)
def fit_func_tan2_fixed_400(P2, b):
    return 1/(-1/400+b*P2)
def fit_func_tan4(P2, a, b, c):
    return 1/(a+b*P2+c*P2*P2)
def fit_func_tan4_fixed_400(P2, b, c):
    return 1/(-1/400+b*P2+c*P2*P2)
def fit_func_tan6(P2, a, b, c, d):
    return 1/(a+b*P2+c*P2*P2+d*P2*P2*P2)

def plot_fits_tan_PS(P_2, fit_result):
    xarr = np.logspace(np.log(min(P_2)[0]*0.9),np.log(max(P_2)[0]*1.1), 500)
    yarr_2 = fit_func_tan2(xarr, fit_result["a2"].median[0], fit_result["b2"].median[0])
    yarr_2p = fit_func_tan2(xarr, fit_result["a2"].median_p[0], fit_result["b2"].median[0])
    yarr_2m = fit_func_tan2(xarr, fit_result["a2"].median_m[0], fit_result["b2"].median[0])
    yarr_2_fixed = fit_func_tan2_fixed_400(xarr, fit_result["b2_fixed"].median[0])
    yarr_2p_fixed = fit_func_tan2_fixed_400(xarr, fit_result["b2_fixed"].median_p[0])
    yarr_2m_fixed = fit_func_tan2_fixed_400(xarr, fit_result["b2_fixed"].median_m[0])
    yarr_4 = fit_func_tan4(xarr, fit_result["a4"].median[0], fit_result["b4"].median[0], fit_result["c4"].median[0])
    yarr_4p = fit_func_tan4(xarr, fit_result["a4"].median_p[0], fit_result["b4"].median[0], fit_result["c4"].median[0])
    yarr_4m = fit_func_tan4(xarr, fit_result["a4"].median_m[0], fit_result["b4"].median[0], fit_result["c4"].median[0])
    yarr_4_fixed = fit_func_tan4_fixed_400(xarr, fit_result["b4_fixed"].median[0], fit_result["c4_fixed"].median[0])
    yarr_4p_fixed = fit_func_tan4_fixed_400(xarr, fit_result["b4_fixed"].median_p[0], fit_result["c4_fixed"].median[0])
    yarr_4m_fixed = fit_func_tan4_fixed_400(xarr, fit_result["b4_fixed"].median_m[0], fit_result["c4_fixed"].median[0])
    yarr_6 = fit_func_tan6(xarr, fit_result["a6"].median[0], fit_result["b6"].median[0], fit_result["c6"].median[0], fit_result["d6"].median[0])
    yarr_6p = fit_func_tan6(xarr, fit_result["a6"].median_p[0], fit_result["b6"].median[0], fit_result["c6"].median[0], fit_result["d6"].median[0])
    yarr_6m = fit_func_tan6(xarr, fit_result["a6"].median_m[0], fit_result["b6"].median[0], fit_result["c6"].median[0], fit_result["d6"].median[0])
    print("a_0(P²): %f + %f - %f"%(fit_result["a_0_2"].median[0],fit_result["a_0_2"].ep[0],fit_result["a_0_2"].em[0]))
    print("a_0(P⁴): %f + %f - %f"%(fit_result["a_0_4"].median[0],fit_result["a_0_4"].ep[0],fit_result["a_0_4"].em[0]))
    print("a_0(P⁶): %f + %f - %f"%(fit_result["a_0_6"].median[0],fit_result["a_0_6"].ep[0],fit_result["a_0_6"].em[0]))
    plt.plot(xarr, yarr_2, color = color_arr[0], label = "$\\mathcal{O}(P^2)$") #: a_0 = %1.3f +  %1.3f -  %1.3f$"%(fit_result["a_0_2"].median[0],fit_result["a_0_2"].ep[0],fit_result["a_0_2"].em[0]))
    plt.fill_between(x=xarr, y1=yarr_2m, y2=yarr_2p, color = color_arr[0], alpha = 0.3)
    # plt.plot(xarr, yarr_2_fixed, color = color_arr[3], label = "$\\mathcal{O}(P^2)$ (fixed)")
    # plt.fill_between(x=xarr, y1=yarr_2m_fixed, y2=yarr_2p_fixed, color = color_arr[3], alpha = 0.3)
    plt.plot(xarr, yarr_4, color = color_arr[1], label = "$\\mathcal{O}(P^4)$") #: a_0 = %1.3f +  %1.3f -  %1.3f$"%(fit_result["a_0_4"].median[0],fit_result["a_0_4"].ep[0],fit_result["a_0_4"].em[0]))
    plt.fill_between(x=xarr, y1=yarr_4m, y2=yarr_4p, color = color_arr[1], alpha = 0.3)
    plt.plot(xarr, yarr_4_fixed, color = color_arr[4], label = "$\\mathcal{O}(P^4)$ (fixed)")
    plt.fill_between(x=xarr, y1=yarr_4m_fixed, y2=yarr_4p_fixed, color = color_arr[4], alpha = 0.3)
    # plt.plot(xarr, yarr_6, color = color_arr[2], label = "$\\mathcal{O}(P^6): a_0 = %1.3f +  %1.3f -  %1.3f$"%(fit_result["a_0_6"].median[0],fit_result["a_0_6"].ep[0],fit_result["a_0_6"].em[0]))
    # plt.fill_between(x=xarr, y1=yarr_6p, y2=yarr_6m, color = color_arr[2], alpha = 0.3)

def plot_P_cot_PS(fit = False, zoom = [0.008,10,-0.11,0.2], title = "", save_app = "", log = ""):
    results = get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot")
    fit_result = get_result_file("phase_shift_fit_P_cot_PS")

    P_cot_PS, P_cot_PS_err, P_2, P_2_err, rho_pi = [[],[],[],[],[]]
    for i in range(len(results)):
        rho_pi.append(results[i]["infos/m_rho_pi"])
        P_2.append(results[i]["P_2_prime"].median)
        P_2_err.append([results[i]["P_2_prime"].ep,results[i]["P_2_prime"].em])
        P_cot_PS.append(results[i]["P_cot_PS_prime"].median[0])
        P_cot_PS_err.append([results[i]["P_cot_PS_prime"].ep,results[i]["P_cot_PS_prime"].em])
    
    norm = matplotlib.colors.Normalize(vmin=min(rho_pi), vmax=max(rho_pi), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    time_color = np.array([(mapper.to_rgba(v)) for v in rho_pi])
    plt.colorbar(mappable=mapper, label = "$\\frac{m_{\\rho}}{m_{\\pi}}$")
    for i in range(len(results)):
        plt.errorbar(x=P_2[i], xerr = P_2_err[i], y=P_cot_PS[i], yerr=P_cot_PS_err[i], color = time_color[i], ls = "", capsize=5, markersize=10)

    if fit:
        plot_fits_P_cot_PS(P_2, fit_result)

    plt.xlim((zoom[0], zoom[1]))
    plt.ylim((zoom[2], zoom[3]))

    plt.axvline(0.4, color = "black", ls="--")
    log_txt = ""
    if log:
        plt.xscale("log")
        log_txt = "_log"
    plt.ylabel("$\\frac{Pcot(\delta)}{m_{\pi}}$")
    plt.xlabel("$\\frac{P^2}{m_{\pi}^2}$")
    plt.legend(fontsize = "x-small")
    plt.title(title)
    plt.grid()
    zoom_txt = ""
    if zoom == plot_tan_PS.__defaults__[0]:
        zoom_txt = "_zoom"
    plt.savefig("plots/P_cot_PS"+log_txt+save_app+zoom_txt+".pdf")
    plt.show()
    plt.clf()

def plot_tan_PS(fit = False, zoom = [0.008,10,-400,100], title = "", save_app = "", log = ""):
    results = get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot")
    fit_result = get_result_file("phase_shift_fit_P_cot_PS")

    tan_PS, tan_PS_err, P_2, P_2_err, rho_pi = [[],[],[],[],[]]
    for i in range(len(results)):
        rho_pi.append(results[i]["infos/m_rho_pi"])
        P_2.append(results[i]["P_2_prime"].median)
        P_2_err.append([results[i]["P_2_prime"].ep,results[i]["P_2_prime"].em])
        tan_PS.append(results[i]["tan_PS/P_prime"].median[0])
        tan_PS_err.append([results[i]["tan_PS/P_prime"].ep,results[i]["tan_PS/P_prime"].em])
    
    norm = matplotlib.colors.Normalize(vmin=min(rho_pi), vmax=max(rho_pi), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    time_color = np.array([(mapper.to_rgba(v)) for v in rho_pi])
    plt.colorbar(mappable=mapper, label = "$\\frac{m_{\\rho}}{m_{\\pi}}$")
    for i in range(len(results)):
        plt.errorbar(x=P_2[i], xerr = P_2_err[i], y=tan_PS[i], yerr=tan_PS_err[i], color = time_color[i], ls = "", capsize=5, markersize=10)

    if fit:
        plot_fits_tan_PS(P_2, fit_result)

    plt.xlim((zoom[0], zoom[1]))
    plt.ylim((zoom[2], zoom[3]))

    plt.title(title)
    plt.axvline(0.4, color = "black", ls="--")
    log_txt = ""
    if log:
        plt.xscale("log")
        log_txt = "_log"
    plt.ylabel("$tan(\delta)\\frac{m_{\pi}}{P}$")
    plt.xlabel("$\\frac{P^2}{m_{\pi}^2}$")
    plt.legend(fontsize = "x-small")
    plt.grid()
    zoom_txt = ""
    if zoom == plot_tan_PS.__defaults__[0]:
        zoom_txt = "_zoom"
    plt.savefig("plots/tan_PS"+log_txt+save_app+zoom_txt+".pdf")
    plt.show()
    plt.clf()







def main():
    # create_all_filenames()
    plot_P_cot_PS(fit = True, title="P cot vs P², Data don't show y-axis crossing", zoom=[-0.04,0.2,-0.06,0.01],save_app="1")
    plot_P_cot_PS(fit = True, title="P cot vs P², Data don't show y-axis crossing", zoom=[-0.05,0.6,-0.15,0.01],save_app="2")
    plot_tan_PS(fit = True, title="Data is consistent with slope getting more and more negative", zoom=[-0.02,1e-1,-500,50])
    plot_tan_PS(fit = True, title="The fitting function demands a constant emerging for small P²\nThis is not happening in the data. Set to -400 here", zoom=[1e-4,1e1,-500,100],log="log")
    plot_P_cot_PS(fit = True, title="In this plot it looks like 1/a_0 is consistent with 0", save_app="fixed", zoom=[1e-4,1e1,-0.15,0.05],log="log")
    # plot_P_cot_PS(fit = True, title="Der Vergleich mit der fixed Methode im cot-plot", save_app="fixed", zoom=[1e-4,1e1,-0.15,0.05])
    # plot_P_cot_PS(fit = True, title="Hier sieht man, dass höhere Ordnungen den Fit nicht verbessern", zoom=[1e-4,1e1,-0.15,0.05],log="log")
    # plot_zeta()
    # plot_dispersion_relation()


def print_phase_shift_data():
    with open("output/phase_shift_data.dat", "w") as file:
        with open("output/phase_shift_data_mathematica.dat", "w") as file_mathematica:
            file.write("P_2\tP_2_err\P_cot\tP_cot_err\ttan/P\ttan/P_err\tm_rho/m_pi\n")
            file_mathematica.write("i\tP_2\tP_2_err\P_cot\tP_cot_err\ttan/P\ttan/P_err\tm_rho/m_pi\n{")
            i = 1
            for result in get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot"):
                rho_pi = result["infos/m_rho_pi"]
                file.write("%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\n"%(result["P_2_prime"].median[0],result["P_2_prime"].e[0],result["P_cot_PS_prime"].median[0],result["P_cot_PS_prime"].e[0],result["tan_PS/P_prime"].median[0],result["P_cot_PS_prime"].e[0],rho_pi))

                rho_pi = result["infos/m_rho_pi"]
                file_mathematica.write("{%i,%.30f,%.30f,%.30f,%.30f,%.30f,%.30f,%.30f},"%(i,result["P_2_prime"].median[0],result["P_2_prime"].e[0],result["P_cot_PS_prime"].median[0],result["P_cot_PS_prime"].e[0],result["tan_PS/P_prime"].median[0],result["P_cot_PS_prime"].e[0],rho_pi))                                  # DONT forget to erase last comma
                i = i+1
            file_mathematica.write("}")

def print_energy_level_data():
    with open("output/energy_level_data.dat", "w") as file:
    # with open("output/energy_level_data_short.dat", "w") as file:
        file.write("gauge_group\tbeta\tm_1\tm_2\tN_T\tN_L\tlevel\tm_rho_pi\tm_pi_inf\tm_pi_inf_err\tE\t\tE_err\n")
        for result in get_results_from_files("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift"):
            for i in range(len(result["E"].median)):
                file.write("%s\t\t"%result["infos/gauge_group"])
                file.write("%1.3f\t"%result["infos/beta"])
                file.write("%1.3f\t"%result["infos/m_1"])
                file.write("%1.3f\t"%result["infos/m_2"])
                file.write("%i\t"%result["N_Ts"].median[i])
                file.write("%i\t"%result["N_Ls"].median[i])
                file.write("%i\t"%result["infos/level"])
                file.write("%e\t"%result["infos/m_rho_pi"])
                file.write("%e\t"%result["infos/m_pi_inf"])
                file.write("%e\t"%result["infos/m_pi_inf_err"])
                file.write("%e\t"%result["E"].median[i])
                file.write("%e\t"%result["E"].e[i])
                file.write("\n")
            file.write("\n")
            
                # file.write("\n___________________________________________________________________________________________________\n")
                
                # print(result.keys(),"\n")
                # file.write("P_2\tP_2_err\P_cot\tP_cot_err\ttan/P\ttan/P_err\tm_rho/m_pi\n")

if __name__ == "__main__":
    # main()
    # print_phase_shift_data()
    print_energy_level_data()












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