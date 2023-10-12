import numpy as np
import error_classes as errcl
import read_HDF5_logfile as HDF_log 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 

color_arr = ["blue", "green", "red", "purple", "lime", "black"]

def infinite_volume(data, args):       
    result = {}
    E_0s = data
    N_Ls = args[0]
    mass_Goldstone = args[1]
    N_Ls_inv = []
    for N_L in N_Ls:
        N_Ls_inv.append(1./N_L)
    def inf_mass_fit_Goldstone(N_L, m_inf, A):
        return m_inf*(1+A*np.exp(-m_inf*N_L)/(m_inf*N_L)**(3/2.))
    def inf_mass_fit_meson(N_L, m_inf, A):
        return m_inf*(1+A*np.exp(-mass_Goldstone*N_L)/(mass_Goldstone*N_L)**(3/2.))
    if mass_Goldstone == 0:
        popt, pcov = curve_fit(f=inf_mass_fit_Goldstone, xdata=N_Ls, ydata=E_0s)
    else:
        popt, pcov = curve_fit(f=inf_mass_fit_meson, xdata=N_Ls, ydata=E_0s)
    result["mass_Goldstone"] = [mass_Goldstone,]
    result["E_0s"] = E_0s
    result["N_Ls"] = N_Ls
    result["N_Ls_inv"] = N_Ls_inv
    result["m_inf"] = [popt[0],]
    result["A_M"] = [popt[1],]
    return result

# def set_errorbar_settings():
#     plt.rcParams["errorbar.capsize"] = 5
#     plt.rcParams["lines.linestyle"] = ""
#     plt.rcParams["lines.markersize"] = 10

# set_errorbar_settings()

################################ CALCULATION ####################################

def main():
    ops = ("pi", "rho", "pipi")
    # ops = ("pi",)
    filelist_list = []
    with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/HDF5_filelist_inifinite_volume", "r") as f:
        for lines in f.readlines():
            if lines[0] != "#":
                filelist_list.append(lines.split())
    for filelist in filelist_list:
        for op in ops:
            E_0s = []
            N_Ls = []
            for files in filelist:
                print(files)
                info = HDF_log.get_info_from_HDF5_logfile(files)
                meas_energlev = errcl.measurement("energy_levels_%s_%s"%(info[7], op), measure_func = None, sampling_args = (None))
                meas_energlev.read_from_HDF()
                # E_0s.append(np.swapaxes(meas_energlev.results["E_0"].sample,0,1)[0])
                E_0s.append(np.swapaxes(meas_energlev.results["E_0_const"].sample,0,1)[0])
                N_Ls.append(info[0])
                if op == "pi":
                    mass_Goldstone = 0
                else:
                    inf_vol = errcl.measurement("infinite_volume_%s_pi"%(info_str), measure_func = None, sampling_args = None)
                    inf_vol.read_from_HDF()
                    mass_Goldstone = inf_vol.results["m_inf"].median[0]
            info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(str(info[6]),info[2],info[3],info[4],info[5])
            inf_vol = errcl.measurement("infinite_volume_%s_%s"%(info_str, op), measure_func = infinite_volume, sampling_args = ("DONT_RESAMPLE",))
            inf_vol.measure(orig_sample=E_0s, args=[N_Ls,mass_Goldstone])
            inf_vol.print_to_HDF()

    for filelist in filelist_list:
        color_ind = 0
        info = HDF_log.get_info_from_HDF5_logfile(filelist[0])
        E_0s = {}
        E_0_eps = {}
        E_0_ems = {}
        m_inf = {}
        m_inf_eps = {}
        m_inf_ems = {}
        A_M = {}
        info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(str(info[6]),info[2],info[3],info[4],info[5])
        for op in ops:
            inf_vol = errcl.measurement("infinite_volume_%s_%s"%(info_str, op), measure_func = None, sampling_args = None)
            inf_vol.read_from_HDF()
            E_0s[op] = inf_vol.results["E_0s"].median
            E_0_eps[op] = inf_vol.results["E_0s"].ep
            E_0_ems[op] = inf_vol.results["E_0s"].em
            m_inf[op] = inf_vol.results["m_inf"].median[0]
            m_inf_eps[op] = inf_vol.results["m_inf"].ep[0]
            m_inf_ems[op] = inf_vol.results["m_inf"].em[0]
            A_M[op] = inf_vol.results["A_M"].median[0]
            mass_Goldstone = inf_vol.results["mass_Goldstone"].median[0]
            N_Ls = inf_vol.results["N_Ls"].median
            N_Ls_inv = inf_vol.results["N_Ls_inv"].median
            def inf_mass_fit_Goldstone(N_L, m_inf, A):
                return m_inf*(1+A*np.exp(-m_inf*N_L)/(m_inf*N_L)**(3/2.))
            def inf_mass_fit_meson(N_L, m_inf, A):
                return m_inf*(1+A*np.exp(-mass_Goldstone*N_L)/(mass_Goldstone*N_L)**(3/2.))
            xarr_inv = np.linspace(1e-10,max(N_Ls_inv), 100)
            xarr = []
            for x in xarr_inv:
                xarr.append(1./x)
            yarr = []
            yarrp = []
            yarrm = []
            for x in xarr:
                if mass_Goldstone == 0:
                    yarr.append(inf_mass_fit_Goldstone(x, m_inf[op], A_M[op]))
                    yarrp.append(inf_mass_fit_Goldstone(x, m_inf[op]+m_inf_eps[op], A_M[op]))
                    yarrm.append(inf_mass_fit_Goldstone(x, m_inf[op]-m_inf_ems[op], A_M[op]))
                else:
                    yarr.append(inf_mass_fit_meson(x, m_inf[op], A_M[op]))
                    yarrp.append(inf_mass_fit_meson(x, m_inf[op]+m_inf_eps[op], A_M[op]))
                    yarrm.append(inf_mass_fit_meson(x, m_inf[op]-m_inf_ems[op], A_M[op]))
            clr = color_arr[color_ind]
            plt.fill_between(x=xarr_inv, y1=yarrp, y2=yarrm, alpha = 0.5, color=clr)
            plt.errorbar(x=N_Ls_inv, y=E_0s[op], yerr=(E_0_eps[op],E_0_ems[op]), label = op, ls = "", capsize=5, markersize=10, color=clr)
            plt.plot(xarr_inv, yarr, c=clr)
            color_ind += 1
        plt.ylabel("E")
        plt.xlabel("1/$N_L$")
        plt.legend()
        plt.grid()
        plt.title(info_str)
        # plt.show()
        plt.savefig("plots/infinite_volume_"+info_str+".pdf")
        plt.clf()



main()





