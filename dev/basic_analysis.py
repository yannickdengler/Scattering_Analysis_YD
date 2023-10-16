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

################################ CALCULATION ####################################

def main():
    filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/HDF5_filelist_full", "str")
    ops = ("pi", "rho", "pipi")
    for filename in filelist:
        info = HDF_log.get_info_from_HDF5_logfile(filename)
        corrs = HDF_log.get_pi_rho_pipi_corr_from_HDF5_logfile(filename)
        for i in range(len(corrs)):
            info["op"] = ops[i]
            basic = errcl.measurement("basic_%s_%s"%(info["info_string"], ops[i]), measure_func = basic_analysis, sampling_args = ("BS_SAMEDIM",1000,1), infos=info)
            basic.measure(orig_sample=np.swapaxes(corrs[i],0,1), args=[None,])
            basic.print_to_HDF()

        # C_plot = {}
        # for i, op in zip(range(3),("pi", "rho", "pipi")):
        #     basic = errcl.measurement("basic_%s_%s"%(info["info_string"], op))
        #     basic.read_from_HDF()

        #     C_plot[op+"_m"] = basic.results["C"].median
        #     C_plot[op+"_ep"] = basic.results["C"].ep
        #     C_plot[op+"_em"] = basic.results["C"].em
        #     plt.errorbar(x=np.arange(len(C_plot[op+"_m"])), y=C_plot[op+"_m"], yerr=(C_plot[op+"_ep"],C_plot[op+"_em"]), label = op)
        # plt.yscale("log")
        # plt.ylabel("C")
        # plt.xlabel("$n_t$")
        # plt.legend()
        # plt.grid()
        # plt.title(info["info_string"])
        # plt.savefig("plots/Corr_"+info["info_string"]+".pdf")
        # plt.clf()

        # m_eff_plot = {}
        # for op in ("pi", "rho", "pipi"):
        #     basic = errcl.measurement("basic_%s_%s"%(info["info_string"], op))
        #     basic.read_from_HDF()
        #     m_eff_plot[op+"_m"] = basic.results["m_eff_impl_deri"].median
        #     m_eff_plot[op+"_ep"] = basic.results["m_eff_impl_deri"].ep
        #     m_eff_plot[op+"_em"] = basic.results["m_eff_impl_deri"].em
        #     plt.errorbar(x=np.arange(len(m_eff_plot[op+"_m"])), y=m_eff_plot[op+"_m"], yerr=(m_eff_plot[op+"_ep"],m_eff_plot[op+"_em"]), label = op)
        # plt.ylabel("C")
        # plt.xlabel("$n_t$")
        # plt.legend()
        # plt.grid()
        # plt.title(info["info_string"])
        # plt.savefig("plots/m_eff_"+info["info_string"]+".pdf")
        # plt.clf()
        
if __name__ == "__main__":
    main()


