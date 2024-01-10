import numpy as np
import error_classes as errcl
import read_HDF5_logfile as HDF_log 
from scipy.optimize import curve_fit
import basic_analysis as basic
import os
import h5py

def calc_const(Correlator,args):                             # Just one Correlator as input
    def corr_fit_fun_cosh(n_t, E_0, A_0, const):
        return A_0*np.cosh((N_T/2.-n_t)*E_0)+const
    result = {}
    fit_limits = args[0]
    N_T = len(Correlator)
    
    xdata = np.arange(N_T)
    popt, pcov = curve_fit(f=corr_fit_fun_cosh, xdata=xdata[fit_limits[0]:fit_limits[1]+1], ydata=Correlator[fit_limits[0]:fit_limits[1]+1])
    result["E_const"] = [abs(popt[0]),]
    result["A_const"] = [popt[1],]
    result["const"] = [popt[2],]
    return result

def calc_E_0_A_0(Correlator,args):                             # Just one Correlator as input
    result = {}
    fit_limits = args[0]
    N_T = len(Correlator)
    C_tilde = np.zeros(len(Correlator)-2)
    for i in range(len(C_tilde)):
        C_tilde[i] = Correlator[i]-Correlator[i+2]
    def corr_fit_func_sinh(n_t, E_0, A_0):
        return A_0*np.sinh((N_T/2.-n_t-1)*E_0)
    
    xdata = np.arange(N_T-1)
    popt, pcov = curve_fit(f=corr_fit_func_sinh, xdata=xdata[fit_limits[0]:fit_limits[1]+1], ydata=C_tilde[fit_limits[0]:fit_limits[1]+1])
    result["E"] = [popt[0],]
    result["A"] = [popt[1],]
    return result

def energy_levels(Correlators, args):                                                      # one correlator, args = size-2 vector of fit_limits
    result = {}
    for key, value in calc_const(Correlators, args).items():
        result[key]=value
    for key, value in calc_E_0_A_0(Correlators, args).items():
        result[key]=value
    for key, value in basic.basic_analysis(Correlators, args).items():
        result[key]=value
    return result

def energy_levels_Fabian(data, args):                                                      # one correlator, args = size-2 vector of fit_limits
    result = {}
    result["E"] = data[:len(data)//2]
    result["A"] = data[len(data)//2:]
    return result

################################ CALCULATION ####################################

def create_all_filenames():
    PATH = "output/HDF5_logfiles/"
    temp = "Scattering_src"
    filelist = os.listdir(PATH)
    resultfile_list = []
    num = len(temp)
    for file in filelist:
        length = len(file)
        if file[:num] == temp:
            resultfile_list.append(file[:length-5])         

    with open("input/filenames_energy_levels_all", "w") as file:
        for filename in resultfile_list:
            file.write(PATH+"%s"%filename+".hdf5\n")

def main():
    filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_energy_levels", "str")
    ops = ("pi", "rho", "pipi")
    for filename in filelist:
        info = HDF_log.get_info_from_HDF5_logfile(filename)
        corrs = HDF_log.get_pi_rho_pipi_corr_from_HDF5_logfile(filename)
        fit_limits = HDF_log.get_fit_limits(filename)
        N_T = info["N_T"]
        for i in range(len(corrs)):
            info["op"] = ops[i]
            energy_lev = errcl.measurement("energy_levels_%s_%s"%(info["info_string"], ops[i]), measure_func = energy_levels, sampling_args = ("BS_SAMEDIM",1000,1), infos=info)
            energy_lev.measure(orig_sample=np.swapaxes(corrs[i],0,1), args=[fit_limits[i],])
            energy_lev.print_to_HDF()

def main_Fabian():
    ops = ("pi", "rho", "pipi")
    PATH = "/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/energy_levels_fabian/"
    filelist = os.listdir(PATH)
    for file in filelist:
        with h5py.File(PATH+file) as file:
            info = HDF_log.get_info_from_Fabian_energy_levels(file)
            for op in ops:
                data = []
                for E in file[op+"/E"][()][:2]:
                    data.append([E,])
                for A in file[op+"/A"][()][:2]:
                    data.append([A,])
                err_arr = []
                for err in file[op+"/Delta_E"][()][:2]:
                    err_arr.append(err)
                for err in file[op+"/Delta_A"][()][:2]:
                    err_arr.append(err)
                energy_lev = errcl.measurement("energy_levels_Fabian_%s_%s"%(info["info_string"], op), measure_func = energy_levels_Fabian, sampling_args = ("GAUSSIAN",err_arr,1000,0), infos=info)
                energy_lev.measure(orig_sample=data, args=None)
                energy_lev.results["E"].e = file[op+"/Delta_E"][()][:2]
                energy_lev.results["A"].e = file[op+"/Delta_A"][()][:2]
                energy_lev.print_to_HDF()
                # energy_lev.print_everything()
                # if op == "pi":
                #     print("Fab: ", "b%1.3fm%1.3fL%iT%i,E1: %1.3f,E2: %1.3f"%(file[op+"/beta"][()],file[op+"/m_1"][()],file[op+"/N_L"][()],file[op+"/N_T"][()],file[op+"/E"][0],file[op+"/E"][1]))
                #     print("Mee: ", "b%1.3fm%1.3fL%iT%i,E1: %1.3f,E2: %1.3f"%(energy_lev.infos["beta"],energy_lev.infos["m_1"],energy_lev.infos["N_L"],energy_lev.infos["N_T"],energy_lev.results["E"].median[0],energy_lev.results["E"].median[1]))
                    

if __name__ == "__main__":
    # create_all_filenames()
    # main()
    main_Fabian()





