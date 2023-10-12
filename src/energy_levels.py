import numpy as np
import error_classes as errcl
import read_HDF5_logfile as HDF_log 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def calc_const(Correlator,args):                             # Just one Correlator as input
    def corr_fit_fun_cosh(n_t, E_0, A_0, const):
        return A_0*np.cosh((N_T/2.-n_t)*E_0)+const
    result = {}
    fit_limits = args[0]
    N_T = len(Correlator)
    
    xdata = np.arange(N_T)
    popt, pcov = curve_fit(f=corr_fit_fun_cosh, xdata=xdata[fit_limits[0]:fit_limits[1]+1], ydata=Correlator[fit_limits[0]:fit_limits[1]+1])
    result["E_0_const"] = [popt[0],]
    result["A_0_const"] = [popt[1],]
    result["const"] = [popt[2],]
    return result

def calc_E_0_A_0(Correlator,args):                             # Just one Correlator as input
    def corr_fit_func_sinh(n_t, E_0, A_0):
        return A_0*np.sinh((N_T/2.-n_t)*E_0)
    result = {}
    fit_limits = args[0]
    N_T = len(Correlator)
    C_tilde = np.zeros(len(Correlator)-1)
    for i in range(len(C_tilde)):
        C_tilde[i] = Correlator[i]-Correlator[i+1]
    
    xdata = np.arange(N_T-1)
    popt, pcov = curve_fit(f=corr_fit_func_sinh, xdata=xdata[fit_limits[0]:fit_limits[1]+1], ydata=C_tilde[fit_limits[0]:fit_limits[1]+1])
    result["E_0"] = [popt[0],]
    result["A_0"] = [popt[1],]
    return result

def energy_levels(Correlators, args):                                                      # one correlator, args = size-2 vector of fit_limits
    result = {}
    # N_T = len(Correlators)
    # xdata = np.arange(N_T)
    for key, value in calc_const(Correlators, args).items():
        result[key]=value
    for key, value in calc_E_0_A_0(Correlators, args).items():
        result[key]=value
    return result

################################ CALCULATION ####################################

def main():
    ops = ("pi", "rho", "pipi")
    filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/HDF5_filelist", "str")

    for filename in filelist:
        print(filename)
        info = HDF_log.get_info_from_HDF5_logfile(filename)
        corrs = HDF_log.get_pi_rho_pipi_corr_from_HDF5_logfile(filename)
        fit_limits = HDF_log.get_fit_limits(filename)
        N_T = info[1]
        def corr_fit_func_sinh(n_t, E_0, A_0):
            return A_0*np.sinh((N_T/2.-n_t)*E_0)
        def corr_fit_fun_cosh(n_t, E_0, A_0, const):
            return A_0*np.cosh((N_T/2.-n_t)*E_0)+const
        for i, corr in zip(range(3), corrs):
            print(ops[i])
            meas_energlev = errcl.measurement("energy_levels_%s_%s"%(info[7], ops[i]), measure_func = energy_levels, sampling_args = ("JK_SAMEDIM",0,0))
            meas_energlev.measure(orig_sample=np.swapaxes(corr,0,1), args=[fit_limits[i],])
            meas_energlev.print_to_HDF()

        C_plot = {}
        for op in ("pi", "rho", "pipi"):
            basic = errcl.measurement("basic_%s_%s"%(info[7], op), measure_func = None, sampling_args = ("JK_SAMEDIM",0,0))
            basic.read_from_HDF()
            meas_energlev = errcl.measurement("energy_levels_%s_%s"%(info[7], op), measure_func = None, sampling_args = ("JK_SAMEDIM",0,0))
            meas_energlev.read_from_HDF()
            C_plot[op+"_m"] = basic.results["C"].median
            C_plot[op+"_ep"] = basic.results["C"].ep
            C_plot[op+"_em"] = basic.results["C"].em
            E_0 = meas_energlev.results["E_0_const"]
            A_0 = meas_energlev.results["A_0_const"]
            const = meas_energlev.results["const"]
            xarr = np.arange(len(C_plot[op+"_m"]))
            plt.errorbar(x=xarr, y=C_plot[op+"_m"], yerr=(C_plot[op+"_ep"],C_plot[op+"_em"]), label = op)
            yarr = corr_fit_fun_cosh(xarr, E_0.median, A_0.median, const.median)
            plt.plot(xarr, yarr)
            yarr_p = corr_fit_fun_cosh(xarr, E_0.median[0] + E_0.ep[0], A_0.median[0], const.median[0])
            yarr_m = corr_fit_fun_cosh(xarr, E_0.median[0] - E_0.em[0], A_0.median[0], const.median[0])
            plt.fill_between(x=xarr, y1=yarr_p, y2=yarr_m, alpha = 0.3)
        plt.yscale("log")
        plt.ylabel("C")
        plt.xlabel("$n_t$")
        plt.legend()
        plt.grid()
        plt.title(info[7])
        plt.savefig("plots/Corr_fit_"+info[7]+".pdf")
        plt.clf()

        m_eff_plot = {}
        for op in ("pi", "rho", "pipi"):
            basic = errcl.measurement("basic_%s_%s"%(info[7], op), measure_func = None, sampling_args = ("JK_SAMEDIM",0,0))
            basic.read_from_HDF()
            meas_energlev = errcl.measurement("energy_levels_%s_%s"%(info[7], op), measure_func = None, sampling_args = ("JK_SAMEDIM",0,0))
            meas_energlev.read_from_HDF()
            m_eff_plot[op+"_m"] = basic.results["m_eff_impl_deri"].median
            m_eff_plot[op+"_ep"] = basic.results["m_eff_impl_deri"].ep
            m_eff_plot[op+"_em"] = basic.results["m_eff_impl_deri"].em
            E_0 = meas_energlev.results["E_0"]
            A_0 = meas_energlev.results["A_0"]
            plt.errorbar(x=np.arange(len(m_eff_plot[op+"_m"])), y=m_eff_plot[op+"_m"], yerr=(m_eff_plot[op+"_ep"],m_eff_plot[op+"_em"]), label = op)
            plt.axhline(y=E_0.median)
            plt.fill_between(x=(0,len(xarr)), y1=(E_0.median[0]+E_0.ep[0],E_0.median[0]+E_0.ep[0]), y2=(E_0.median[0]-E_0.em[0],E_0.median[0]-E_0.em[0]))
        plt.ylabel("C")
        plt.xlabel("$n_t$")
        plt.legend()
        plt.grid()
        plt.title(info[7])
        plt.savefig("plots/m_eff_fit_"+info[7]+".pdf")
        plt.clf()


main()





