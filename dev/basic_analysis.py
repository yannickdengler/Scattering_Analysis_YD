import numpy as np
from scipy.optimize import bisect
import error_classes as errcl
import read_HDF5_logfile as HDF_log 
import matplotlib.pyplot as plt


def calc_corr(Correlators, args):                               # Correlators [Ops][N_T] 
    result = {}      
    result["C"] = Correlators                              
    return result

def calc_corr_tilde(Correlators, args):                               # Correlators [Ops][N_T] 
    result = {}               
    result["C_tilde"] = []
    for i in range(len(Correlators)-2):                                        # 2 oder 1, I have to decide
        result["C_tilde"].append(Correlators[i]-Correlators[i+2])
    return result

def calc_eff_mass_log(Correlators, args):                               # Correlators [Ops][N_T] 
    result = {}              
    result["m_eff_log"] = []
    for i in range(len(Correlators)-1):
        m_eff = np.log(Correlators[i])/np.log(Correlators[i+1])
        if np.isnan(m_eff) or np.isinf(m_eff):
            result["m_eff_log"].append(0)
        else:
            result["m_eff_log"].append(m_eff)
    return result

def calc_eff_mass_impl(Correlators, args):                               # Correlators [Ops][N_T] 
    result = {}                                    
    def zero_eff_mass(eff_mass, ratio, index):
        return np.cosh(eff_mass*(T_2-index))/np.cosh(eff_mass*(T_2-(index+1))) - ratio                                               # {results}[result_array]
    
    result["m_eff_impl"] = []
    for i in range(len(Correlators)-1):                                        # 2 oder 1, I have to decide
        ratio = Correlators[i]/Correlators[i+1]
        T_2 = (len(Correlators)-2)//2
        result["m_eff_impl"].append(bisect(f=zero_eff_mass, a=1e-30, b=100, args = (ratio,i)))
    return result

def calc_eff_mass_impl_deri(Correlators, args):                               # Correlators [Ops][N_T] 
    result = {}                                    
    def zero_eff_mass(eff_mass, ratio, index):
        return np.sinh(eff_mass*(T_2-index))/np.sinh(eff_mass*(T_2-(index+1))) - ratio                                               # {results}[result_array]
    result["m_eff_impl_deri"] = []
    for i in range(len(Correlators)-3):                                        # 2 oder 1, I have to decide
        ratio = (Correlators[i]-Correlators[i+2])/(Correlators[i+1]-Correlators[i+3])
        T_2 = (len(Correlators)-2)//2
        if (T_2-i) == 0 or (T_2-(i+1)) == 0:                    # Only happens for sinh. For cosh both values are well-defined
            # result["m_eff_impl_deri_"+Op].append(float("inf"))
            result["m_eff_impl_deri"].append(0)
        else:
            res = bisect(f=zero_eff_mass, a=1e-30, b=1000, args = (ratio,i))
            if np.isnan(res):
                result["m_eff_impl_deri"].append(0)
            else:
                result["m_eff_impl_deri"].append(bisect(f=zero_eff_mass, a=1e-30, b=1000, args = (ratio,i)))
    return result
    
def calc_convexity(Correlators, args):                          # Correlators [Ops][N_T] 
    result = {}                                                 # {results}[result_array]
    result["convexity"] = []
    for i in range(len(Correlators)-2):
        result["convexity"].append((Correlators[i]-2*Correlators[i+1]+Correlators[i+2])/4)
    return result

def basic_analysis(Correlators, args):                                                      # Give C_pi, C_rho, C_pipi
    result = {}
    for key, value in calc_corr(Correlators, None).items():
        result[key]=value
    for key, value in calc_corr_tilde(Correlators, None).items():
        result[key]=value
    for key, value in calc_eff_mass_log(Correlators, None).items():
        result[key]=value
    for key, value in calc_eff_mass_impl_deri(Correlators, None).items():
        result[key]=value
    for key, value in calc_convexity(Correlators, None).items():
        result[key]=value
    return result

# def set_errorbar_settings():
#     plt.rcParams["errorbar.capsize"] = 5
#     plt.rcParams["lines.linestyle"] = ""
#     plt.rcParams["lines.markersize"] = 10

# set_errorbar_settings()


################################ CALCULATION ####################################

def main():
    filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/HDF5_filelist_full", "str")

    for filename in filelist:
        info = HDF_log.get_info_from_HDF5_logfile(filename)
        corrs = HDF_log.get_pi_rho_pipi_corr_from_HDF5_logfile(filename)
        for op, corr in zip(("pi", "rho", "pipi"), corrs):
            basic = errcl.measurement("basic_%s_%s"%(info[7], op), measure_func = basic_analysis, sampling_args = ("BS_SAMEDIM",600,1))
            basic.measure(orig_sample=np.swapaxes(corr,0,1), args=[None,])
            basic.print_to_HDF()

        C_plot = {}
        for i, op in zip(range(3),("pi", "rho", "pipi")):
            basic = errcl.measurement("basic_%s_%s"%(info[7], op), measure_func = basic_analysis, sampling_args = None)
            basic.read_from_HDF()

            C_plot[op+"_m"] = basic.results["C"].median
            C_plot[op+"_ep"] = basic.results["C"].ep
            C_plot[op+"_em"] = basic.results["C"].em
            plt.errorbar(x=np.arange(len(C_plot[op+"_m"])), y=C_plot[op+"_m"], yerr=(C_plot[op+"_ep"],C_plot[op+"_em"]), label = op)
        plt.yscale("log")
        plt.ylabel("C")
        plt.xlabel("$n_t$")
        plt.legend()
        plt.grid()
        plt.title(info[7])
        plt.savefig("plots/Corr_"+info[7]+".pdf")
        plt.clf()

        m_eff_plot = {}
        for op in ("pi", "rho", "pipi"):
            basic = errcl.measurement("basic_%s_%s"%(info[7], op), measure_func = None, sampling_args = None)
            basic.read_from_HDF()
            m_eff_plot[op+"_m"] = basic.results["m_eff_impl_deri"].median
            m_eff_plot[op+"_ep"] = basic.results["m_eff_impl_deri"].ep
            m_eff_plot[op+"_em"] = basic.results["m_eff_impl_deri"].em
            plt.errorbar(x=np.arange(len(m_eff_plot[op+"_m"])), y=m_eff_plot[op+"_m"], yerr=(m_eff_plot[op+"_ep"],m_eff_plot[op+"_em"]), label = op)
        plt.ylabel("C")
        plt.xlabel("$n_t$")
        plt.legend()
        plt.grid()
        plt.title(info[7])
        plt.savefig("plots/m_eff_"+info[7]+".pdf")
        plt.clf()

main()





        # print(len(corrs),len(corrs[0]),len(corrs[0][0]))
        # corr_pi = corrs[0]
        # corr_5 = np.swapaxes(corr_pi,0,1)[5]
        # print(len(corr_5))
        # num = len(corr_5)
        # mean_np = np.mean(corr_5)
        # mean_std = np.std(corr_5)
        # var = 0
        # for x in corr_5:
        #     var += (x-mean_np)**2/(num-1)
        # std = np.sqrt(var)
        # std_err = std/np.sqrt(num)

        # # print(mean_std)
        # # print(std)

        # # print(mean_std/np.sqrt(num-1))
        # # print(std_err)


        # # exit()

        # counter = 0
        # for x in corr_5:
        #     if x < mean_np+std_err and x > mean_np-std_err:
        #         counter += 1
        # print(counter, num, counter/num)

        # print(num)
        # print("std:", mean_np, std, std/mean_np)
        # print("std_err:", mean_np, std_err, std_err/mean_np)
        # basic = errcl.measurement("basic_%s_%s_2000BS"%(info[7], op), measure_func = basic_analysis, sampling_args = None)
        # basic.read_from_HDF()
        # print("68%:", basic.results["C"].median[5], basic.results["C"].e[5], basic.results["C"].e[5]/basic.results["C"].median[5])
        # print("BS:", basic.results["C"].median[5], basic.results["C"].e_BS[5], basic.results["C"].e_BS[5]/basic.results["C"].median[5])
        # basic = errcl.measurement("basic_%s_%s_500BS"%(info[7], op), measure_func = basic_analysis, sampling_args = None)
        # basic.read_from_HDF()
        # print("68%:", basic.results["C"].median[5], basic.results["C"].e[5], basic.results["C"].e[5]/basic.results["C"].median[5])
        # print("BS:", basic.results["C"].median[5], basic.results["C"].e_BS[5], basic.results["C"].e_BS[5]/basic.results["C"].median[5])


        # basic_JK = errcl.measurement("basic_%s_%s"%(info[7], op), measure_func = basic_analysis, sampling_args = None)
        # basic_JK.read_from_HDF()
        # print("JK:", basic_JK.results["C"].mean[5], basic_JK.results["C"].e_JK[5], basic_JK.results["C"].e_JK[5]/basic_JK.results["C"].mean[5])

