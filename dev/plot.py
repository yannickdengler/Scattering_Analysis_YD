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
        # print(resultfile)
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
                # plt.show()
                plt.savefig("plots/m_eff_fit_Fabian_2_"+infos["info_string"]+".pdf")
                plt.clf()

def inf_mass_fit_Goldstone(N_L, m_inf, A):
    return m_inf*(1+A*np.exp(-m_inf*N_L)/(m_inf*N_L)**(3/2.))
def inf_mass_fit_meson(N_L, m_inf, A, mass_Goldstone):
    return m_inf*(1+A*np.exp(-mass_Goldstone*N_L)/(mass_Goldstone*N_L)**(3/2.))

def plot_infinite_volume():
    ops = ("pi", "rho", "pipi")
    # ops = ("pi",)
    # resultfile_list = get_result_files("infinite_volume")
    resultfile_list = get_result_files("infinite_volume_Fabian_level_0")
    for resultfile in resultfile_list:
        color_ind = 0
        if resultfile[len(resultfile)-3:] == "_pi":
            file = resultfile[:len(resultfile)-2]
            # inf_vol = errcl.measurement(file+"pi")
            # inf_vol.read_from_HDF()
            # infos = inf_vol.infos
            # print(infos["level"])
            for op in ops:
                inf_vol = errcl.measurement(file+op)
                inf_vol.read_from_HDF()
                infos = inf_vol.infos
                E_0 = inf_vol.results["E"].mean
                E_0_e = inf_vol.results["E"].e
                # E_0 = inf_vol.results["E"].median
                # E_0_ep = inf_vol.results["E"].ep
                # E_0_em = inf_vol.results["E"].em
                m_inf = inf_vol.results["m_inf"].mean[0]
                m_inf_e = inf_vol.results["m_inf"].e[0]
                # m_inf = inf_vol.results["m_inf"].median[0]
                # m_inf_ep = inf_vol.results["m_inf"].ep[0]
                # m_inf_em = inf_vol.results["m_inf"].em[0]
                A_M = inf_vol.results["A_M"].mean[0]
                # A_M = inf_vol.results["A_M"].median[0]
                mass_Goldstone = inf_vol.results["mass_Goldstone"].mean[0]
                # mass_Goldstone = inf_vol.results["mass_Goldstone"].median[0]
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
                        yarrp.append(inf_mass_fit_Goldstone(x, m_inf+m_inf_e, A_M))
                        yarrm.append(inf_mass_fit_Goldstone(x, m_infmm_inf_e, A_M))
                        # yarrp.append(inf_mass_fit_Goldstone(x, m_inf+m_inf_ep, A_M))
                        # yarrm.append(inf_mass_fit_Goldstone(x, m_inf-m_inf_em, A_M))
                    else:
                        yarr.append(inf_mass_fit_meson(x, m_inf, A_M, mass_Goldstone))
                        yarrp.append(inf_mass_fit_meson(x, m_inf+m_inf_e, A_M, mass_Goldstone))
                        yarrm.append(inf_mass_fit_meson(x, m_inf-m_inf_e, A_M, mass_Goldstone))
                        # yarrp.append(inf_mass_fit_meson(x, m_inf+m_inf_ep, A_M, mass_Goldstone))
                        # yarrm.append(inf_mass_fit_meson(x, m_inf-m_inf_em, A_M, mass_Goldstone))
                clr = color_arr[color_ind]
                plt.fill_between(x=xarr_inv, y1=yarrp, y2=yarrm, alpha = 0.5, color=clr)
                plt.errorbar(x=N_Ls_inv, y=E_0, yerr=E_0_e, label = op, ls = "", capsize=5, markersize=10, color=clr)
                plt.plot(xarr_inv, yarr, c=clr)
                color_ind += 1
            plt.ylabel("E")
            plt.xlabel("1/$N_L$")
            plt.legend()
            plt.grid()
            plt.title(infos["info_string"])
            # plt.show()
            plt.savefig("plots/infinite_volume_"+infos["info_string"]+".pdf")
            plt.clf()

def compare_energy_levels():
    x_axis_counter = 0
    ops = ("pi", "rho", "pipi")
    resultfile_list = get_result_files("infinite_volume_Fabian")
    print(resultfile_list)
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



def main():
    # plot_basic_analysis()
    # plot_energy_levels()
    plot_energy_levels()
    # plot_infinite_volume()
    # compare_energy_levels()
    # plot_phase_shift()


if __name__ == "__main__":
    main()