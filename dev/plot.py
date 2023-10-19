import os
import matplotlib.pyplot as plt
import read_HDF5_logfile as HDF_log
import error_classes as errcl 
import numpy as np
import matplotlib
from matplotlib import cm



color_arr = ["blue", "green", "red", "purple", "lime", "black", "grey", "orange", "olive", "skyblue", "fuchsia", "peru", "firebrick"]

def get_result_files(temp):
    PATH = "output/result_files/"
    filelist = os.listdir(PATH)
    resultfile_list = []
    num = len(temp)
    for file in filelist:
        length = len(file)
        if file[:num] == temp:
            resultfile_list.append(file[:length-5])                                 # -5 removes the ".hdf5"
    return resultfile_list

def plot_basic_analysis(plot_corr = True, plot_m_eff = True):
    ops = ("pi", "rho", "pipi")
    resultfile_list = get_result_files("basic_Scattering")
    infos = {}
    for resultfile in resultfile_list:
        if resultfile[len(resultfile)-3:] == "_pi":
            file = resultfile[:len(resultfile)-2]
            if plot_corr:
                C = {}
                C_err = {}
                for op in ops:
                    basic = errcl.measurement(file+op)
                    basic.read_from_HDF()
                    C[op] = basic.results["C"].median
                    C_err[op] = (basic.results["C"].ep,basic.results["C"].em)
                    plt.errorbar(x=np.arange(len(C[op])), y=C[op], yerr=C_err[op], label = op)
                    infos = basic.infos
                plt.yscale("log")
                plt.ylabel("C")
                plt.xlabel("$n_t$")
                plt.legend()
                plt.grid()
                plt.title(infos["info_string"])
                # plt.show()
                plt.savefig("plots/Corr_"+infos["info_string"]+".pdf")
                plt.clf()
            if plot_m_eff:
                m_eff = {}
                m_eff_err = {}
                for op in ops:
                    basic = errcl.measurement(file+op)
                    basic.read_from_HDF()
                    m_eff[op] = basic.results["m_eff_impl_deri"].median
                    m_eff_err[op] = (basic.results["m_eff_impl_deri"].ep,basic.results["m_eff_impl_deri"].em)
                    plt.errorbar(x=np.arange(len(m_eff[op])), y=m_eff[op], yerr=m_eff_err[op], label = op)
                    infos = basic.infos
                plt.yscale("log")
                plt.ylabel("m_eff")
                plt.xlabel("$n_t$")
                plt.legend()
                plt.grid()
                plt.title(infos["info_string"])
                # plt.show()
                plt.savefig("plots/m_eff_"+infos["info_string"]+".pdf")
                plt.clf()

def corr_fit_func_sinh(n_t, E_0, A_0, N_T):
    return A_0*np.sinh((N_T/2.-n_t-1)*E_0)
def corr_fit_fun_cosh(n_t, E_0, A_0, const, N_T):
    return A_0*np.cosh((N_T/2.-n_t)*E_0)+const

def plot_energy_levels(plot_corr = True, plot_m_eff = True):
    ops = ("pi", "rho", "pipi")
    resultfile_list = get_result_files("energy_levels")
    for resultfile in resultfile_list:
        if resultfile[len(resultfile)-3:] == "_pi":
            file = resultfile[:len(resultfile)-2]
            if plot_corr:
                for op in ops:
                    meas_energlev = errcl.measurement(file+op)
                    meas_energlev.read_from_HDF()
                    infos = meas_energlev.infos
                    C = meas_energlev.results["C"].median
                    C_err = (meas_energlev.results["C"].ep,meas_energlev.results["C"].em)
                    E_0 = meas_energlev.results["E_const"]
                    A_0 = meas_energlev.results["A_const"]
                    const = meas_energlev.results["const"]
                    xarr = np.linspace(0, len(C),100)
                    yarr = corr_fit_fun_cosh(xarr, E_0.median[0], A_0.median[0], const.median[0], infos["N_T"])
                    yarr_p = corr_fit_fun_cosh(xarr, E_0.median[0] + E_0.ep[0], A_0.median[0], const.median[0], infos["N_T"])
                    yarr_m = corr_fit_fun_cosh(xarr, E_0.median[0] - E_0.em[0], A_0.median[0], const.median[0], infos["N_T"])
                    cor = color_arr[ops.index(op)]
                    plt.errorbar(x=np.arange(len(C)), y=C, yerr=C_err, label = op, color = cor)
                    plt.plot(xarr, yarr, color = cor)
                    plt.fill_between(x=xarr, y1=yarr_p, y2=yarr_m, alpha = 0.3, color = cor)
                plt.yscale("log")
                plt.ylabel("C")
                plt.xlabel("$n_t$")
                plt.legend()
                plt.grid()
                plt.title(infos["info_string"])
                # plt.show()
                plt.savefig("plots/Corr_fit_"+infos["info_string"]+".pdf")
                plt.clf()
            if plot_m_eff:
                for op in ops:
                    meas_energlev = errcl.measurement(file+op)
                    meas_energlev.read_from_HDF()
                    infos = meas_energlev.infos
                    m_eff = meas_energlev.results["m_eff_impl_deri"].median
                    m_eff_err = (meas_energlev.results["m_eff_impl_deri"].ep,meas_energlev.results["m_eff_impl_deri"].em)
                    E_0 = meas_energlev.results["E"]
                    A_0 = meas_energlev.results["A"]
                    xarr = np.linspace(0, len(m_eff),100)
                    cor = color_arr[ops.index(op)]
                    plt.errorbar(x=np.arange(len(m_eff)), y=m_eff, yerr=m_eff_err, label = op, color = cor)
                    plt.axhline(y=E_0.median[0], color = cor)
                    plt.fill_between(x=(0,len(m_eff)), y1=(E_0.median[0]+E_0.ep[0],E_0.median[0]+E_0.ep[0]), y2=(E_0.median[0]-E_0.em[0],E_0.median[0]-E_0.em[0]), color = cor)
                plt.xlim((0,infos["N_T"]))
                plt.ylabel("m_eff")
                plt.xlabel("$n_t$")
                plt.legend()
                plt.grid()
                plt.title(infos["info_string"])
                # plt.show()
                plt.savefig("plots/m_eff_fit_"+infos["info_string"]+".pdf")
                plt.clf()

def corr_fit_fun_cosh2(n_t, E_0, A_0, E_1, A_1, N_T):
    return A_0*np.cosh((N_T/2.-n_t)*E_0)+ A_1*np.cosh((N_T/2.-n_t)*E_1)

def plot_energy_levels_Fabian(plot_corr = True, plot_m_eff = True):
    ops = ("pi", "rho", "pipi")
    resultfile_list = get_result_files("energy_levels_Fabian")
    for resultfile in resultfile_list:
        print(resultfile)
        if resultfile[len(resultfile)-3:] == "_pi":
            file = resultfile[:len(resultfile)-2]
            if plot_corr:
                for op in ops:
                    meas_energlev = errcl.measurement(file+op)
                    meas_energlev.read_from_HDF()
                    infos = meas_energlev.infos
                    basic = errcl.measurement("basic_"+infos["info_string"]+"_%s"%op)
                    if basic.file_exists():
                        basic.read_from_HDF()
                        C = basic.results["C"].median
                        C_err = (basic.results["C"].ep,basic.results["C"].em)
                        E = meas_energlev.results["E"]
                        A = meas_energlev.results["A"]
                        # const = meas_energlev.results["const"]
                        xarr = np.linspace(0, len(C),100)
                        yarr = corr_fit_fun_cosh2(xarr, E.median[0], A.median[0], E.median[1], A.median[1], infos["N_T"])
                        yarr_p = corr_fit_fun_cosh2(xarr, E.median[0] + E.ep[0], A.median[0], E.median[1] + E.ep[1], A.median[1], infos["N_T"])
                        yarr_m = corr_fit_fun_cosh2(xarr, E.median[0] - E.em[0], A.median[0], E.median[1] - E.em[1], A.median[1], infos["N_T"])
                        cor = color_arr[ops.index(op)]
                        plt.errorbar(x=np.arange(len(C)), y=C, yerr=C_err, label = op, color = cor)
                        plt.plot(xarr, yarr, color = cor)
                        plt.fill_between(x=xarr, y1=yarr_p, y2=yarr_m, alpha = 0.3, color = cor)
                plt.yscale("log")
                plt.ylabel("C")
                plt.xlabel("$n_t$")
                plt.legend()
                plt.grid()
                plt.title(infos["info_string"])
                # plt.show()
                plt.savefig("plots/Corr_fit_Fabian_2_"+infos["info_string"]+".pdf")
                plt.clf()
            if plot_m_eff:
                for op in ops:
                    meas_energlev = errcl.measurement(file+op)
                    meas_energlev.read_from_HDF()
                    infos = meas_energlev.infos
                    basic = errcl.measurement("basic_"+infos["info_string"]+"_%s"%op)
                    if basic.file_exists():
                        basic.read_from_HDF()
                        m_eff = basic.results["m_eff_impl_deri"].median
                        m_eff_err = (basic.results["m_eff_impl_deri"].ep,basic.results["m_eff_impl_deri"].em)
                        E_0 = meas_energlev.results["E"]
                        A_0 = meas_energlev.results["A"]
                        xarr = np.linspace(0, len(m_eff),100)
                        cor = color_arr[ops.index(op)]
                        plt.errorbar(x=np.arange(len(m_eff)), y=m_eff, yerr=m_eff_err, label = op, color = cor)
                        plt.axhline(y=E_0.median[0], color = cor)
                        plt.fill_between(x=(0,len(m_eff)), y1=(E_0.median[0]+E_0.ep[0],E_0.median[0]+E_0.ep[0]), y2=(E_0.median[0]-E_0.em[0],E_0.median[0]-E_0.em[0]), color = cor)
                plt.xlim((0,infos["N_T"]))
                plt.ylabel("m_eff")
                plt.xlabel("$n_t$")
                plt.legend()
                plt.grid()
                plt.title(infos["info_string"])
                plt.show()
                plt.savefig("plots/m_eff_fit_Fabian_2_"+infos["info_string"]+".pdf")
                plt.clf()

def inf_mass_fit_Goldstone(N_L, m_inf, A):
    return m_inf*(1+A*np.exp(-m_inf*N_L)/(m_inf*N_L)**(3/2.))
def inf_mass_fit_meson(N_L, m_inf, A, mass_Goldstone):
    return m_inf*(1+A*np.exp(-mass_Goldstone*N_L)/(mass_Goldstone*N_L)**(3/2.))

def plot_infinite_volume():
    ops = ("pi", "rho", "pipi")
    # resultfile_list = get_result_files("infinite_volume")
    resultfile_list = get_result_files("infinite_volume_Fabian")
    for resultfile in resultfile_list:
        color_ind = 0
        if resultfile[len(resultfile)-3:] == "_pi":
            file = resultfile[:len(resultfile)-2]
            for op in ops:
                inf_vol = errcl.measurement(file+op)
                inf_vol.read_from_HDF()
                infos = inf_vol.infos
                E_0 = inf_vol.results["E"].median
                E_0_ep = inf_vol.results["E"].ep
                E_0_em = inf_vol.results["E"].em
                m_inf = inf_vol.results["m_inf"].median[0]
                m_inf_ep = inf_vol.results["m_inf"].ep[0]
                m_inf_em = inf_vol.results["m_inf"].em[0]
                A_M = inf_vol.results["A_M"].median[0]
                mass_Goldstone = inf_vol.results["mass_Goldstone"].median[0]
                N_Ls = inf_vol.results["N_Ls"].median
                N_Ls_inv = inf_vol.results["N_Ls_inv"].median
                xarr_inv = np.linspace(1e-10,max(N_Ls_inv), 100)
                xarr = []
                for x in xarr_inv:
                    xarr.append(1./x)
                yarr = []
                yarrp = []
                yarrm = []
                for x in xarr:
                    if mass_Goldstone == 0:
                        yarr.append(inf_mass_fit_Goldstone(x, m_inf, A_M))
                        yarrp.append(inf_mass_fit_Goldstone(x, m_inf+m_inf_ep, A_M))
                        yarrm.append(inf_mass_fit_Goldstone(x, m_inf-m_inf_em, A_M))
                    else:
                        yarr.append(inf_mass_fit_meson(x, m_inf, A_M, mass_Goldstone))
                        yarrp.append(inf_mass_fit_meson(x, m_inf+m_inf_ep, A_M, mass_Goldstone))
                        yarrm.append(inf_mass_fit_meson(x, m_inf-m_inf_em, A_M, mass_Goldstone))
                clr = color_arr[color_ind]
                plt.fill_between(x=xarr_inv, y1=yarrp, y2=yarrm, alpha = 0.5, color=clr)
                plt.errorbar(x=N_Ls_inv, y=E_0, yerr=(E_0_ep,E_0_em), label = op, ls = "", capsize=5, markersize=10, color=clr)
                plt.plot(xarr_inv, yarr, c=clr)
                color_ind += 1
            plt.ylabel("E")
            plt.xlabel("1/$N_L$")
            plt.legend()
            plt.grid()
            plt.title(infos["info_string"])
            plt.show()
            plt.savefig("plots/infinite_volume_"+infos["info_string"]+".pdf")
            plt.clf()

def compare_energy_levels():
    x_axis_counter = 0
    ops = ("pi", "rho", "pipi")
    resultfile_list = get_result_files("infinite_volume_Scattering")
    for resultfile in resultfile_list:
        if resultfile[len(resultfile)-3:] == "_pi":
            file = resultfile[:len(resultfile)-2]
            for op in ops:
                temp = file[16:]
                print(temp)
                inf_vol = errcl.measurement("infinite_volume_"+temp+op)
                inf_vol.read_from_HDF()
                inf_vol_const = errcl.measurement("infinite_volume_const_"+temp+op)
                inf_vol_const.read_from_HDF()
                inf_vol_Fabian = errcl.measurement("infinite_volume_Fabian_"+temp+op)
                inf_vol_Fabian.read_from_HDF()
                # inf_vol_Fabian2 = errcl.measurement("infinite_volume_Fabian_level2_"+temp+op)
                # inf_vol_Fabian2.read_from_HDF()
                infos = inf_vol.infos
                N_Ls = inf_vol.results["N_Ls"].median
                N_Ls_Fabian = inf_vol_Fabian.results["N_Ls"].median
                for N_L in N_Ls:
                    plt.errorbar(x=x_axis_counter, y=inf_vol.results["E"].median[N_Ls.index(N_L)], yerr = inf_vol.results["E"].e[N_Ls.index(N_L)], color = "red", elinewidth=10)
                    x_axis_counter += 1
                    plt.errorbar(x=x_axis_counter, y=inf_vol_const.results["E"].median[N_Ls.index(N_L)], yerr = inf_vol_const.results["E"].e[N_Ls.index(N_L)], color = "green", elinewidth=10)
                    x_axis_counter += 1
                    plt.errorbar(x=x_axis_counter, y=inf_vol_Fabian.results["E"].median[N_Ls_Fabian.index(N_L)], yerr = inf_vol_Fabian.results["E"].e[N_Ls_Fabian.index(N_L)], color = "blue", elinewidth=10)
                    x_axis_counter += 2
                x_axis_counter += 1
                plt.errorbar(x=x_axis_counter, y=inf_vol.results["m_inf"].median[0], yerr = inf_vol.results["m_inf"].e[0], color = "darkred", elinewidth=10)
                x_axis_counter += 1
                plt.errorbar(x=x_axis_counter, y=inf_vol_const.results["m_inf"].median[0], yerr = inf_vol_const.results["m_inf"].e[0], color = "darkgreen", elinewidth=10)
                x_axis_counter += 1
                plt.errorbar(x=x_axis_counter, y=inf_vol_Fabian.results["m_inf"].median[0], yerr = inf_vol_Fabian.results["m_inf"].e[0], color = "darkblue", elinewidth=10)
   
                x_axis_counter += 2.5
                plt.axvline(x_axis_counter)
                x_axis_counter += 2.5
            plt.title("%1.3f %1.3f"%(infos["beta"],infos["m_1"]))
            plt.show()

def plot_phase_shift():
    ensemble_list = []
    color_ind = -1


    y = []
    yerr = []
    x = []
    xerr = []
    rho_pi = []

    # resultfile_list = get_result_files("phase_shift")
    # resultfile_list = get_result_files("phase_shift_const")
    resultfile_list = get_result_files("phase_shift_Fabian")
    for resultfile in resultfile_list:
        phase_shift = errcl.measurement(resultfile)
        phase_shift.read_from_HDF()
        info = phase_shift.infos
        pi_energy_lev = errcl.measurement("energy_levels_Fabian_Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i_pi"%(info["isospin_channel"],info["gauge_group"],info["beta"],info["m_1"],info["m_2"],info["N_T"],info["N_L"]))
        pi_energy_lev.read_from_HDF()
        rho_energy_lev = errcl.measurement("energy_levels_Fabian_Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i_rho"%(info["isospin_channel"],info["gauge_group"],info["beta"],info["m_1"],info["m_2"],info["N_T"],info["N_L"]))
        rho_energy_lev.read_from_HDF()
        pi_energy_lev.print_everything()
        rho_pi.append(rho_energy_lev.results["E"].median[0]/pi_energy_lev.results["E"].median[0])
        infos = phase_shift.infos
        N_T = infos["N_T"]
        N_L = infos["N_L"]
        info_str = "Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(infos["isospin_channel"],infos["gauge_group"],infos["beta"],infos["m_1"],infos["m_2"])
        ensemble_str = "beta=%1.3f, m1/2=%1.3f"%(infos["beta"],infos["m_1"])
        y.append(phase_shift.results["P_cot_PS_prime"].median)
        yerr.append([phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em])
        x.append(phase_shift.results["P_2_prime"].median)
        xerr.append([phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em])
        # if not (ensemble_str in ensemble_list):
        #     color_ind += 1
        #     ensemble_list.append(ensemble_str)
        #     plt.errorbar(x=phase_shift.results["P_2_prime"].median, xerr = (phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em), y=phase_shift.results["P_cot_PS_prime"].median, yerr=(phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em), label = ensemble_str, ls = "", capsize=5, markersize=10, color = color_arr[color_ind])
        # else:
        #     plt.errorbar(x=phase_shift.results["P_2_prime"].median, xerr = (phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em), y=phase_shift.results["P_cot_PS_prime"].median, yerr=(phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em), ls = "", capsize=5, markersize=10, color = color_arr[color_ind])


    norm = matplotlib.colors.Normalize(vmin=min(rho_pi), vmax=max(rho_pi), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    time_color = np.array([(mapper.to_rgba(v)) for v in rho_pi])
    for i in range(len(yerr)):
        plt.errorbar(x=x[i], xerr = xerr[i], y=y[i], yerr=yerr[i], color = time_color[i], ls = "", capsize=5, markersize=10)
        
    plt.xscale("log")
    plt.ylabel("$\\frac{Pcot(\delta)}{m_{\pi}}$")
    plt.xlabel("$\\frac{P^2}{m_{\pi}}$")
    plt.legend(fontsize = "xx-small")
    plt.xlim((0.008, 10))
    plt.ylim((-0.11,0.2))
    plt.grid()
    # plt.title(info[7])
    plt.savefig("plots/P_cot_PS_.pdf")
    plt.show()
    plt.clf()

def fit_func_p2(P2, a, b):
    return a + b*P2
def fit_func_p4(P2, a, b, c):
    return a + b*P2 + c*P2*P2
def fit_func_p6(P2, a, b, c, d):
    return a + b*P2 + c*P2*P2 + d*P2*P2*P2

def plot_phase_shift_fit(zoom = [0.008,10,-0.11,0.2]):
    ensemble_list = []
    color_ind = -1


    y = []
    yerr = []
    x = []
    xerr = []
    rho_pi = []

    # resultfile_list = get_result_files("phase_shift")
    # resultfile_list = get_result_files("phase_shift_const")
    resultfile_list = get_result_files("phase_shift_Fabian")
    for resultfile in resultfile_list:
        phase_shift = errcl.measurement(resultfile)
        phase_shift.read_from_HDF()
        info = phase_shift.infos
        pi_energy_lev = errcl.measurement("energy_levels_Fabian_Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i_pi"%(info["isospin_channel"],info["gauge_group"],info["beta"],info["m_1"],info["m_2"],info["N_T"],info["N_L"]))
        pi_energy_lev.read_from_HDF()
        rho_energy_lev = errcl.measurement("energy_levels_Fabian_Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i_rho"%(info["isospin_channel"],info["gauge_group"],info["beta"],info["m_1"],info["m_2"],info["N_T"],info["N_L"]))
        rho_energy_lev.read_from_HDF()
        # pi_energy_lev.print_everything()
        rho_pi.append(rho_energy_lev.results["E"].median[0]/pi_energy_lev.results["E"].median[0])
        infos = phase_shift.infos
        N_T = infos["N_T"]
        N_L = infos["N_L"]
        info_str = "Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(infos["isospin_channel"],infos["gauge_group"],infos["beta"],infos["m_1"],infos["m_2"])
        ensemble_str = "beta=%1.3f, m1/2=%1.3f"%(infos["beta"],infos["m_1"])
        y.append(phase_shift.results["P_cot_PS_prime"].median)
        yerr.append([phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em])
        x.append(phase_shift.results["P_2_prime"].median)
        xerr.append([phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em])
        # if info["level"] == 0:
        #     y.append(phase_shift.results["P_cot_PS_prime"].median)
        #     yerr.append([phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em])
        #     x.append(phase_shift.results["P_2_prime"].median)
        #     xerr.append([phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em])
        # if info["level"] == 1:
        #     y1.append(phase_shift.results["P_cot_PS_prime"].median)
        #     y1err.append([phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em])
        #     x1.append(phase_shift.results["P_2_prime"].median)
        #     x1err.append([phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em])



    phase_shift_fit = errcl.measurement("phase_shift_fit")
    phase_shift_fit.read_from_HDF()
    xarr = np.logspace(np.log(min(x)[0]*0.9),np.log(max(x)[0]*1.1), 500)
    # xarr = np.linspace(min(min(x))*0.9, max(max(x))*1.1, 500)
    yarr_2 = fit_func_p2(xarr, phase_shift_fit.results["a2"].median[0], phase_shift_fit.results["b2"].median[0])
    yarr_2p = fit_func_p2(xarr, phase_shift_fit.results["a2"].median_p[0], phase_shift_fit.results["b2"].median[0])
    yarr_2m = fit_func_p2(xarr, phase_shift_fit.results["a2"].median_m[0], phase_shift_fit.results["b2"].median[0])
    yarr_4 = fit_func_p4(xarr, phase_shift_fit.results["a4"].median[0], phase_shift_fit.results["b4"].median[0], phase_shift_fit.results["c4"].median[0])
    yarr_4p = fit_func_p4(xarr, phase_shift_fit.results["a4"].median_p[0], phase_shift_fit.results["b4"].median[0], phase_shift_fit.results["c4"].median[0])
    yarr_4m = fit_func_p4(xarr, phase_shift_fit.results["a4"].median_m[0], phase_shift_fit.results["b4"].median[0], phase_shift_fit.results["c4"].median[0])
    yarr_6 = fit_func_p6(xarr, phase_shift_fit.results["a6"].median[0], phase_shift_fit.results["b6"].median[0], phase_shift_fit.results["c6"].median[0], phase_shift_fit.results["d6"].median[0])
    yarr_6p = fit_func_p6(xarr, phase_shift_fit.results["a6"].median_p[0], phase_shift_fit.results["b6"].median[0], phase_shift_fit.results["c6"].median[0], phase_shift_fit.results["d6"].median[0])
    yarr_6m = fit_func_p6(xarr, phase_shift_fit.results["a6"].median_m[0], phase_shift_fit.results["b6"].median[0], phase_shift_fit.results["c6"].median[0], phase_shift_fit.results["d6"].median[0])

    print("a_0(P²): %f + %f - %f"%(phase_shift_fit.results["a_0_2"].median[0],phase_shift_fit.results["a_0_2"].ep[0],phase_shift_fit.results["a_0_2"].em[0]))
    print("a_0(P⁴): %f + %f - %f"%(phase_shift_fit.results["a_0_4"].median[0],phase_shift_fit.results["a_0_4"].ep[0],phase_shift_fit.results["a_0_4"].em[0]))
    print("a_0(P⁶): %f + %f - %f"%(phase_shift_fit.results["a_0_6"].median[0],phase_shift_fit.results["a_0_6"].ep[0],phase_shift_fit.results["a_0_6"].em[0]))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.plot(xarr, yarr_2, color = color_arr[0], label = "$\\mathcal{O}(P^2): a_0 = %1.3f +  %1.3f -  %1.3f$"%(phase_shift_fit.results["a_0_2"].median[0],phase_shift_fit.results["a_0_2"].ep[0],phase_shift_fit.results["a_0_2"].em[0]))
    plt.fill_between(x=xarr, y1=yarr_2m, y2=yarr_2p, color = color_arr[0], alpha = 0.3)
    plt.plot(xarr, yarr_4, color = color_arr[1], label = "$\\mathcal{O}(P^4): a_0 = %1.3f +  %1.3f -  %1.3f$"%(phase_shift_fit.results["a_0_4"].median[0],phase_shift_fit.results["a_0_4"].ep[0],phase_shift_fit.results["a_0_4"].em[0]))
    plt.fill_between(x=xarr, y1=yarr_4m, y2=yarr_4p, color = color_arr[1], alpha = 0.3)
    plt.plot(xarr, yarr_6, color = color_arr[2], label = "$\\mathcal{O}(P^6): a_0 = %1.3f +  %1.3f -  %1.3f$"%(phase_shift_fit.results["a_0_6"].median[0],phase_shift_fit.results["a_0_6"].ep[0],phase_shift_fit.results["a_0_6"].em[0]))
    plt.fill_between(x=xarr, y1=yarr_6p, y2=yarr_6m, color = color_arr[2], alpha = 0.3)

    norm = matplotlib.colors.Normalize(vmin=min(rho_pi), vmax=max(rho_pi), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
    time_color = np.array([(mapper.to_rgba(v)) for v in rho_pi])
    plt.colorbar(mappable=mapper, label = "$\\frac{m_{\\rho}}{m_{\\pi}}$")
    for i in range(len(yerr)):
        plt.errorbar(x=x[i], xerr = xerr[i], y=y[i], yerr=yerr[i], color = time_color[i], ls = "", capsize=5, markersize=10)
    plt.xlim((zoom[0], zoom[1]))
    plt.ylim((zoom[2], zoom[3]))

    plt.axvline(0.4, color = "black", ls="--")
    plt.xscale("log")
    plt.ylabel("$\\frac{Pcot(\delta)}{m_{\pi}}$")
    plt.xlabel("$\\frac{P^2}{m_{\pi}^2}$")
    plt.legend(fontsize = "x-small")
    plt.grid()
    # plt.title(info[7])
    if zoom == plot_phase_shift_fit.__defaults__[0]:
        plt.savefig("plots/P_cot_PS_fit.pdf")
    else:
        plt.savefig("plots/P_cot_PS_fit_zoom.pdf")
    plt.show()
    plt.clf()

def print_phase_shift_data():
    with open("output/phase_shift_data.dat", "w") as file:
        file.write("P_cot_PS\tP_cot_PS_err\tP_2\tP_2_err\tm_rho/m_pi\n")

        # resultfile_list = get_result_files("phase_shift")
        # resultfile_list = get_result_files("phase_shift_const")
        resultfile_list = get_result_files("phase_shift_Fabian")
        for resultfile in resultfile_list:
            phase_shift = errcl.measurement(resultfile)
            phase_shift.read_from_HDF()
            info = phase_shift.infos
            pi_energy_lev = errcl.measurement("energy_levels_Fabian_Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i_pi"%(info["isospin_channel"],info["gauge_group"],info["beta"],info["m_1"],info["m_2"],info["N_T"],info["N_L"]))
            pi_energy_lev.read_from_HDF()
            rho_energy_lev = errcl.measurement("energy_levels_Fabian_Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i_rho"%(info["isospin_channel"],info["gauge_group"],info["beta"],info["m_1"],info["m_2"],info["N_T"],info["N_L"]))
            rho_energy_lev.read_from_HDF()
            rho_pi = rho_energy_lev.results["E"].median[0]/pi_energy_lev.results["E"].median[0]
            file.write("%e\t%e\t%e\t%e\t%e\n"%(phase_shift.results["P_cot_PS_prime"].median[0],phase_shift.results["P_cot_PS_prime"].e[0],phase_shift.results["P_2_prime"].median[0],phase_shift.results["P_2_prime"].e[0],rho_pi))





def main():
    # plot_basic_analysis()
    # plot_energy_levels()
    # plot_energy_levels_Fabian()
    # plot_infinite_volume()
    # compare_energy_levels()
    # plot_phase_shift()
    plot_phase_shift_fit()
    plot_phase_shift_fit(zoom=[0.01,0.6,-0.11, 0.01])
    # print_phase_shift_data()


if __name__ == "__main__":
    main()