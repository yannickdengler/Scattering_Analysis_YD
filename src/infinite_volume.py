import numpy as np
import error_classes as errcl
import read_HDF5_logfile as HDF_log 
from scipy.optimize import curve_fit 
import os

color_arr = ["blue", "green", "red", "purple", "lime", "black"]

def infinite_volume(data, args):       
    result = {}
    E_0s = data
    N_Ls = args[0]
    mass_Goldstone_input = args[1]
    N_Ts = args[2]
    mass_Goldstone = 0
    N_Ls_inv = []
    for N_L in N_Ls:
        N_Ls_inv.append(1./N_L)
    def inf_mass_fit_Goldstone(N_L, m_inf, A):
        return m_inf*(1+A*np.exp(-m_inf*N_L)/(m_inf*N_L)**(3/2.))
    def inf_mass_fit_meson(N_L, m_inf, A):
        return m_inf*(1+A*np.exp(-mass_Goldstone_input*N_L)/(mass_Goldstone_input*N_L)**(3/2.))
    if mass_Goldstone_input == 0:
        popt, pcov = curve_fit(f=inf_mass_fit_Goldstone, xdata=N_Ls, ydata=E_0s)
        mass_Goldstone = popt[0]
    else:
        popt, pcov = curve_fit(f=inf_mass_fit_meson, xdata=N_Ls, ydata=E_0s)
        mass_Goldstone = mass_Goldstone_input
    result["mass_Goldstone_input"] = [mass_Goldstone_input,]
    result["mass_Goldstone"] = [mass_Goldstone,]
    result["E"] = E_0s
    result["N_Ts"] = N_Ts
    result["N_Ls"] = N_Ls
    result["N_Ls_inv"] = N_Ls_inv
    result["m_inf"] = [popt[0],]
    result["A_M"] = [popt[1],]
    return result
    
################################ CALCULATION ####################################

def create_all_filenames():
    PATH = "output/result_files/"
    tmp = "energy_levels_Fabian"
    filelist = os.listdir(PATH)
    resultfile_list = []
    num = len(tmp)
    for file in filelist:
        length = len(file)
        if file[:num] == tmp:
            resultfile_list.append(file[:length-5])         

    beta_m_str = ""
    beta_m_str_list = []
    filelist_list = []
    for file in resultfile_list:
        if file[len(file)-5:] == "_pipi":
            meas_energlev = errcl.measurement(file)                     # without hdf5
            meas_energlev.read_from_HDF()
            info = meas_energlev.infos
            beta_m_str = "beta%1.3f_m1%1.3f_m2%1.3f"%(info["beta"],info["m_1"],info["m_2"])
            if not (beta_m_str in beta_m_str_list):
                beta_m_str_list.append(beta_m_str)
                filelist_list.append([])
            filelist_list[beta_m_str_list.index(beta_m_str)].append(file[:len(file)-5])

    with open("input/filenames_infinite_volume_all", "w") as filestream:
        for filelist in filelist_list:
            for file in filelist:
                filestream.write("%s\t"%file)
            filestream.write("\n")

def main():
    ops = ("pi", "rho", "pipi")
    filelist_list = []
    with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_infinite_volume", "r") as f:
        for lines in f.readlines():
            if lines[0] != "#":
                filelist_list.append(lines.split())
    for filelist in filelist_list:
        for op in ops:
            E_0 = []
            E_0_const = []
            N_Ls = []
            N_Ts = []
            for files in filelist:
                info = HDF_log.get_info_from_HDF5_logfile(files)
                info["op"] = op
                meas_energlev = errcl.measurement("energy_levels_%s_%s"%(info["info_string"], op))
                meas_energlev.read_from_HDF()
                E_0_const.append(np.swapaxes(meas_energlev.results["E_0_const"].sample,0,1)[0])
                E_0.append(np.swapaxes(meas_energlev.results["E_0"].sample,0,1)[0])
                N_Ls.append(info["N_L"])
                N_Ts.append(info["N_T"])
                if op == "pi":
                    mass_Goldstone = 0
                else:
                    inf_vol = errcl.measurement("infinite_volume_%s_pi"%(info_str))
                    inf_vol.read_from_HDF()
                    mass_Goldstone = inf_vol.results["m_inf"].median[0]
            info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%("I"+str(info["isospin_channel"]),info["gauge_group"],info["beta"],info["m_1"],info["m_2"])
            inf_vol_const = errcl.measurement("infinite_volume_const_%s_%s"%(info_str, op), measure_func = infinite_volume, sampling_args = ("DONT_RESAMPLE",), infos=info)
            inf_vol_const.measure(orig_sample=E_0_const, args=[N_Ls,mass_Goldstone,N_Ts])
            inf_vol_const.print_to_HDF()
            inf_vol = errcl.measurement("infinite_volume_%s_%s"%(info_str, op), measure_func = infinite_volume, sampling_args = ("DONT_RESAMPLE",), infos=info)
            inf_vol.measure(orig_sample=E_0, args=[N_Ls,mass_Goldstone,N_Ts])
            inf_vol.print_to_HDF()

def main_Fabian():
    ops = ("pi", "rho", "pipi")
    filelist_list = []
    with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_infinite_volume_all", "r") as f:
        for lines in f.readlines():
            if lines[0] != "#":
                filelist_list.append(lines.split())

    for filelist in filelist_list:
        if len(filelist) > 1:
            for op in ops:
                E_0 = []
                E_1 = []
                N_Ls = []
                N_Ts = []
                for files in filelist:
                    print("files: ", files)
                    meas_energlev = errcl.measurement(files+"_"+op)
                    meas_energlev.read_from_HDF()
                    info = meas_energlev.infos                
                    info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%("I"+str(info["isospin_channel"]),info["gauge_group"],info["beta"],info["m_1"],info["m_2"])
                    info["op"] = op
                    info["info_str"] = info_str
                    E_0.append(np.swapaxes(meas_energlev.results["E"].sample,0,1)[0])
                    E_1.append(np.swapaxes(meas_energlev.results["E"].sample,0,1)[1])
                    N_Ls.append(info["N_L"])
                    N_Ts.append(info["N_T"])
                    if op == "pi":
                        mass_Goldstone = 0
                    else:
                        inf_vol = errcl.measurement("infinite_volume_Fabian_level_0_%s_pi"%(info_str))
                        inf_vol.read_from_HDF()
                        mass_Goldstone = inf_vol.results["m_inf"].median[0]
                inf_vol1 = errcl.measurement("infinite_volume_Fabian_level_0_%s_%s"%(info_str, op), measure_func = infinite_volume, sampling_args = ("DONT_RESAMPLE",), infos=info)
                info["level"] = 0
                inf_vol1.measure(orig_sample=E_0, args=[N_Ls,mass_Goldstone,N_Ts])
                inf_vol1.print_to_HDF()
                inf_vol2 = errcl.measurement("infinite_volume_Fabian_level_1_%s_%s"%(info_str, op), measure_func = infinite_volume, sampling_args = ("DONT_RESAMPLE",), infos=info)
                info["level"] = 1
                inf_vol2.measure(orig_sample=E_1, args=[N_Ls,mass_Goldstone,N_Ts])
                inf_vol2.print_to_HDF()



if __name__ == "__main__":
    # create_all_filenames()
    # main()
    main_Fabian()



