import numpy as np
import error_classes as errcl
import read_HDF5_logfile as HDF_log 
import matplotlib.pyplot as plt
# from scipy.optimize import curve_fit 
import generalizedzeta as gz

import h5py
import os

color_arr = ["blue", "green", "red", "purple", "lime", "black", "grey", "orange", "olive", "skyblue", "fuchsia", "peru", "firebrick"]


def generalized_mom(E_pipi, m_meson, L):
    if E_pipi < 2*m_meson:
        return 0
    return 2*np.arcsin(np.sqrt(0.5*(np.cosh(E_pipi/2.)-np.cosh(m_meson))))

def calc_phase_shift(data, args):       
    result = {}
    E_pipi = data[0]
    N_L = args[0]
    mass_Goldstone = args[1]
    N_L_inv = 1./N_L
    P = generalized_mom(E_pipi, mass_Goldstone, N_L)
    q = P*N_L/(2*np.pi)
    Zeta = gz.Zeta(q**2)
    tan_PS = np.pi**(3/2.)*q*Zeta
    # print(mass_Goldstone, E_pipi,P,q,Zeta, tan_PS)
    # if np.isinf(tan_PS) or np.isnan(tan_PS):
    #     tan_PS = 0
    M = 16*np.pi*tan_PS/np.sqrt(tan_PS**2+1)
    result["N_L"] = np.asarray((N_L,))
    result["N_L_inv"] = np.asarray((N_L_inv,))
    result["mass_Goldstone"] = np.asarray((mass_Goldstone,))
    result["E_pipi"] = np.asarray((E_pipi,))
    result["E_pipi_prime"] = np.asarray((E_pipi/mass_Goldstone,))
    result["s_pipi"] = np.asarray((E_pipi**2,))
    result["s_pipi_prime"] = np.asarray(((E_pipi/mass_Goldstone)**2,))
    result["P"] = np.asarray((P,))
    result["P_prime"] = np.asarray((P/mass_Goldstone,))
    result["P_2"] = np.asarray((P**2,))
    result["P_2_prime"] = np.asarray(((P/mass_Goldstone)**2,))
    result["q"] = np.asarray((q,))
    result["q_2"] = np.asarray((q**2,))
    result["tan_PS"] = np.asarray((tan_PS,))
    result["PS"] = np.asarray((np.arctan(tan_PS),))
    result["cot_PS"] = np.asarray((1./tan_PS,))
    result["P_cot_PS"] = np.asarray((P/tan_PS,))
    result["P_cot_PS_prime"] = np.asarray((P/(tan_PS*mass_Goldstone),))
    result["P3_cot_PS_s"] = np.asarray((P**3/(E_pipi*tan_PS*mass_Goldstone**2),))
    result["P3_cot_PS_s_prime"] = np.asarray((P/(E_pipi*tan_PS*mass_Goldstone**2),))
    result["Zeta"] = np.asarray((Zeta,))
    result["M"] = np.asarray((M,))
    result["M_2"] = np.asarray((M*M,))
    return result


################################ CALCULATION ####################################

def main():


    ############################ FABIANS DATA ################

    # PATH = "/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/energy_levels_fabian/"
    # filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/HDF5_filelist_phase_shft", "str")
    # filelist_Fabian = os.listdir(PATH)
    # for filename, filename_Fabian in zip(filelist, filelist_Fabian):
    #     print(filename)
    #     info = HDF_log.get_info_from_HDF5_logfile(filename)
    #     N_T = info[1]
    #     N_L = info[0]
    #     info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(str(info[6]),info[2],info[3],info[4],info[5])
    #     inf_vol = errcl.measurement("infinite_volume_Fabian_%s_pi"%(info_str), measure_func = None, sampling_args = None)
    #     inf_vol.read_from_HDF()
    #     mass_Goldstone = inf_vol.results["m_inf"].median[0]
    #     inf_vol_pipi = errcl.measurement("infinite_volume_Fabian_%s_pipi"%(info_str), measure_func = None, sampling_args = None)
    #     inf_vol_pipi.read_from_HDF()
    #     E_0s = [np.swapaxes(inf_vol_pipi.results["E_0s"].sample,0,1)[inf_vol_pipi.results["N_Ls"].median.index(N_L)],]
    #     if (np.mean(E_0s)*0.99 > 2*mass_Goldstone) and (np.mean(E_0s)*1.01 < 4*mass_Goldstone):
    #         phase_shift = errcl.measurement("phase_shift_Fabian_%s"%(info[7]), measure_func = calc_phase_shift, sampling_args = ["DONT_RESAMPLE",])
    #         phase_shift.measure(orig_sample=E_0s, args=(N_L,mass_Goldstone))
    #         phase_shift.print_to_HDF()
    #     inf_vol_pipi = errcl.measurement("infinite_volume_Fabian_level2_%s_pipi"%(info_str), measure_func = None, sampling_args = None)
    #     inf_vol_pipi.read_from_HDF()
    #     E_0s = [np.swapaxes(inf_vol_pipi.results["E_0s"].sample,0,1)[inf_vol_pipi.results["N_Ls"].median.index(N_L)],]
    #     if (np.mean(E_0s)*0.99 > 2*mass_Goldstone) and (np.mean(E_0s)*1.01 < 4*mass_Goldstone):
    #         phase_shift = errcl.measurement("phase_shift_Fabian_level2_%s"%(info[7]), measure_func = calc_phase_shift, sampling_args = ["DONT_RESAMPLE",])
    #         phase_shift.measure(orig_sample=E_0s, args=(N_L,mass_Goldstone))
    #         phase_shift.print_to_HDF()


    # PATH = "/home/dengler_yannick/Documents/Scattering_Analysis_YD/output/result_files/"
    # filelist = os.listdir(PATH)

    # result_listYD = []
    # result_list1 = []
    # result_list2 = []
    # for file in filelist:
    #     # print(file[:39])
    #     templ = "phase_shift_Scattering_I2_SP(4)_"
    #     if file[:32] == templ:
    #         result_listYD.append(file[:len(file)-5])
    #     templ = "phase_shift_Fabian_Scattering_I2_SP(4)_"
    #     if file[:39] == templ:
    #         result_list1.append(file[:len(file)-5])
    #     templ = "phase_shift_Fabian_level2_Scattering_I2_SP(4)_"
    #     if file[:46] == templ:
    #         result_list2.append(file[:len(file)-5])

    # for file in result_listYD:
    #     print(file)
    #     if file == result_listYD[0]:
    #         plt.scatter(1,1,color = color_arr[0],label ="Yannick")
    #     phase_shift = errcl.measurement(file, measure_func = calc_phase_shift, sampling_args = None)
    #     phase_shift.read_from_HDF()
    #     plt.errorbar(x=phase_shift.results["P_prime"].median, xerr = (phase_shift.results["P_prime"].ep,phase_shift.results["P_prime"].em), y=phase_shift.results["tan_PS"].median, yerr=(phase_shift.results["tan_PS"].ep,phase_shift.results["tan_PS"].em), ls = "", capsize=5, markersize=10, color = color_arr[0])
    # for file in result_list1:
    #     if file == result_list1[0]:
    #         plt.scatter(1,1,color = color_arr[1],label ="Fabian level 1")
    #     phase_shift = errcl.measurement(file, measure_func = calc_phase_shift, sampling_args = None)
    #     phase_shift.read_from_HDF()
    #     plt.errorbar(x=phase_shift.results["P_prime"].median, xerr = (phase_shift.results["P_prime"].ep,phase_shift.results["P_prime"].em), y=phase_shift.results["tan_PS"].median, yerr=(phase_shift.results["tan_PS"].ep,phase_shift.results["tan_PS"].em), ls = "", capsize=5, markersize=10, color = color_arr[1])
    # for file in result_list2:
    #     if file == result_list2[0]:
    #         plt.scatter(1,1,color = color_arr[2],label ="Fabian level 2")
    #     phase_shift = errcl.measurement(file, measure_func = calc_phase_shift, sampling_args = None)
    #     phase_shift.read_from_HDF()
    #     plt.errorbar(x=phase_shift.results["P_prime"].median, xerr = (phase_shift.results["P_prime"].ep,phase_shift.results["P_prime"].em), y=phase_shift.results["tan_PS"].median, yerr=(phase_shift.results["tan_PS"].ep,phase_shift.results["tan_PS"].em), ls = "", capsize=5, markersize=10, color = color_arr[2])
    
    
    
    # plt.xscale("log")
    # plt.ylabel("tan($\delta$)")
    # plt.xlabel("$\\frac{P}{m_{\pi}}$")
    # plt.legend(fontsize = "xx-small")
    # plt.grid()
    # # plt.savefig("plots/tan_PS_Fabian"+".pdf")
    # plt.savefig("plots/tan_PS_Fabian_level2"+".pdf")
    # plt.show()
    # plt.clf()




    # for file in result_listYD:
    #     if file == result_listYD[0]:
    #         plt.scatter(1,1,color = color_arr[0],label ="Yannick")
    #     phase_shift = errcl.measurement(file, measure_func = calc_phase_shift, sampling_args = None)
    #     phase_shift.read_from_HDF()
    #     plt.errorbar(x=phase_shift.results["P_2_prime"].median, xerr = (phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em), y=phase_shift.results["P_cot_PS_prime"].median, yerr=(phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em), ls = "", capsize=5, markersize=10, color = color_arr[0])
    # for file in result_list1:
    #     if file == result_list1[0]:
    #         plt.scatter(1,1,color = color_arr[1],label ="Fabian level 1")
    #     phase_shift = errcl.measurement(file, measure_func = calc_phase_shift, sampling_args = None)
    #     phase_shift.read_from_HDF()
    #     plt.errorbar(x=phase_shift.results["P_2_prime"].median, xerr = (phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em), y=phase_shift.results["P_cot_PS_prime"].median, yerr=(phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em), ls = "", capsize=5, markersize=10, color = color_arr[1])
    # for file in result_list2:
    #     if file == result_list2[0]:
    #         plt.scatter(1,1,color = color_arr[2],label ="Fabian level 2")
    #     phase_shift = errcl.measurement(file, measure_func = calc_phase_shift, sampling_args = None)
    #     phase_shift.read_from_HDF()
    #     plt.errorbar(x=phase_shift.results["P_2_prime"].median, xerr = (phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em), y=phase_shift.results["P_cot_PS_prime"].median, yerr=(phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em), ls = "", capsize=5, markersize=10, color = color_arr[2])
    # plt.xscale("log")
    # plt.ylabel("$\\frac{Pcot(\delta)}{m_{\pi}}$")
    # plt.xlabel("$\\frac{P^2}{m_{\pi}}$")
    # plt.legend(fontsize = "xx-small")
    # plt.grid()
    # # plt.savefig("plots/P_cot_PS_Fabian"+".pdf")
    # plt.savefig("plots/P_cot_PS_Fabian_level2"+".pdf")
    # plt.show()
    # plt.clf()




    # for file in result_listYD:
    #     if file == result_listYD[0]:
    #         plt.scatter(1,1,color = color_arr[0],label ="Yannick")
    #     phase_shift = errcl.measurement(file, measure_func = calc_phase_shift, sampling_args = None)
    #     phase_shift.read_from_HDF()
    #     plt.errorbar(x=phase_shift.results["s_pipi_prime"].median, xerr = (phase_shift.results["s_pipi_prime"].ep,phase_shift.results["s_pipi_prime"].em), y=phase_shift.results["PS"].median, yerr=(phase_shift.results["PS"].ep,phase_shift.results["PS"].em), ls = "", capsize=5, markersize=10, color = color_arr[0])
    # for file in result_list1:    
    #     if file == result_list1[0]:    
    #         plt.scatter(1,1,color = color_arr[1],label ="Fabian level 1")
    #     phase_shift = errcl.measurement(file, measure_func = calc_phase_shift, sampling_args = None)
    #     phase_shift.read_from_HDF()
    #     plt.errorbar(x=phase_shift.results["s_pipi_prime"].median, xerr = (phase_shift.results["s_pipi_prime"].ep,phase_shift.results["s_pipi_prime"].em), y=phase_shift.results["PS"].median, yerr=(phase_shift.results["PS"].ep,phase_shift.results["PS"].em), ls = "", capsize=5, markersize=10, color = color_arr[1])
    # for file in result_list2:
    #     if file == result_list2[0]:
    #         plt.scatter(1,1,color = color_arr[2],label ="Fabian level 2")
    #     phase_shift = errcl.measurement(file, measure_func = calc_phase_shift, sampling_args = None)
    #     phase_shift.read_from_HDF()
    #     plt.errorbar(x=phase_shift.results["s_pipi_prime"].median, xerr = (phase_shift.results["s_pipi_prime"].ep,phase_shift.results["s_pipi_prime"].em), y=phase_shift.results["PS"].median, yerr=(phase_shift.results["PS"].ep,phase_shift.results["PS"].em), ls = "", capsize=5, markersize=10, color = color_arr[2])
    # # plt.xscale("log")
    # plt.ylabel("$\delta$")
    # plt.xlabel("$\\frac{s}{m_{\pi}^2}$")
    # plt.legend(fontsize = "xx-small")
    # plt.grid()
    # # plt.savefig("plots/PS_Fabian"+".pdf")
    # plt.savefig("plots/PS_Fabian_level2"+".pdf")
    # plt.show()
    # plt.clf()




    # for file in result_listYD:
    #     if file == result_listYD[0]:
    #         plt.scatter(1,1,color = color_arr[0],label ="Yannick")
    #     phase_shift = errcl.measurement(file, measure_func = calc_phase_shift, sampling_args = None)
    #     phase_shift.read_from_HDF()
    #     plt.errorbar(x=phase_shift.results["s_pipi_prime"].median, xerr = (phase_shift.results["s_pipi_prime"].ep,phase_shift.results["s_pipi_prime"].em), y=phase_shift.results["P3_cot_PS_s_prime"].median, yerr=(phase_shift.results["P3_cot_PS_s_prime"].ep,phase_shift.results["P3_cot_PS_s_prime"].em), ls = "", capsize=5, markersize=10, color = color_arr[0])
    # for file in result_list1:
    #     if file == result_list1[0]:
    #         plt.scatter(1,1,color = color_arr[1],label ="Fabian level 1")
    #     phase_shift = errcl.measurement(file, measure_func = calc_phase_shift, sampling_args = None)
    #     phase_shift.read_from_HDF()
    #     plt.errorbar(x=phase_shift.results["s_pipi_prime"].median, xerr = (phase_shift.results["s_pipi_prime"].ep,phase_shift.results["s_pipi_prime"].em), y=phase_shift.results["P3_cot_PS_s_prime"].median, yerr=(phase_shift.results["P3_cot_PS_s_prime"].ep,phase_shift.results["P3_cot_PS_s_prime"].em), ls = "", capsize=5, markersize=10, color = color_arr[1])
    # for file in result_list2:
    #     if file == result_list2[0]:
    #         plt.scatter(1,1,color = color_arr[2],label ="Fabian level 2")
    #     phase_shift = errcl.measurement(file, measure_func = calc_phase_shift, sampling_args = None)
    #     phase_shift.read_from_HDF()
    #     plt.errorbar(x=phase_shift.results["s_pipi_prime"].median, xerr = (phase_shift.results["s_pipi_prime"].ep,phase_shift.results["s_pipi_prime"].em), y=phase_shift.results["P3_cot_PS_s_prime"].median, yerr=(phase_shift.results["P3_cot_PS_s_prime"].ep,phase_shift.results["P3_cot_PS_s_prime"].em), ls = "", capsize=5, markersize=10, color = color_arr[2])
    # # plt.xscale("log")
    # plt.ylabel("$\\frac{P^3cot(\delta)}{s^{0.5}m_{\pi}}$")
    # plt.xlabel("$\\frac{s}{m_{\pi}^2}$")
    # plt.legend(fontsize = "xx-small")
    # plt.grid()
    # # plt.savefig("plots/PS_Fabian"+".pdf")
    # plt.savefig("plots/resonance_Fabian_level2"+".pdf")
    # plt.show()
    # plt.clf()

    ############################ FABIANS DATA END ################

    ############################ MY DATA ################


    filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/HDF5_filelist_phase_shift", "str")
    for filename in filelist:
        print(filename)
        info = HDF_log.get_info_from_HDF5_logfile(filename)
        N_T = info["N_T"]
        N_L = info["N_L"]
        info_str = "Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(info["isospin_channel"],info["gauge_group"],info["beta"],info["m_1"],info["m_2"])
        inf_vol = errcl.measurement("infinite_volume_%s_pi"%(info_str))
        # inf_vol = errcl.measurement("infinite_volume_const_%s_pi"%(info_str))
        inf_vol.read_from_HDF()
        mass_Goldstone = inf_vol.results["m_inf"].median[0]
        energy_levels = errcl.measurement("energy_levels_%s_pipi"%(info["info_string"]))
        energy_levels.read_from_HDF()
        E_0s = np.swapaxes(energy_levels.results["E_0"].sample,0,1)
        if energy_levels.results["E_0"].median[0]*0.99 > 2*mass_Goldstone:
        # E_0s = np.swapaxes(energy_levels.results["E_0_const"].sample,0,1)
        # if energy_levels.results["E_0_const"].median[0]*0.99 > 2*mass_Goldstone:
            phase_shift = errcl.measurement("phase_shift_%s"%(info["info_string"]), measure_func = calc_phase_shift, sampling_args = ["DONT_RESAMPLE",])
            # phase_shift = errcl.measurement("phase_shift_const_%s"%(info["info_string"]), measure_func = calc_phase_shift, sampling_args = ["DONT_RESAMPLE",])
            phase_shift.measure(orig_sample=E_0s, args=(N_L,mass_Goldstone))
            phase_shift.print_to_HDF()

    for filename in filelist:
        print(filename)
        info = HDF_log.get_info_from_HDF5_logfile(filename)
        N_T = info["N_T"]
        N_L = info["N_L"]
        info_str = "Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(info["isospin_channel"],info["gauge_group"],info["beta"],info["m_1"],info["m_2"])
        # inf_vol = errcl.measurement("infinite_volume_%s_pi"%(info_str))
        inf_vol = errcl.measurement("infinite_volume_const_%s_pi"%(info_str))
        inf_vol.read_from_HDF()
        mass_Goldstone = inf_vol.results["m_inf"].median[0]
        energy_levels = errcl.measurement("energy_levels_%s_pipi"%(info["info_string"]))
        energy_levels.read_from_HDF()
        # E_0s = np.swapaxes(energy_levels.results["E_0"].sample,0,1)
        # if energy_levels.results["E_0"].median[0]*0.99 > 2*mass_Goldstone:
        E_0s = np.swapaxes(energy_levels.results["E_0_const"].sample,0,1)
        if energy_levels.results["E_0_const"].median[0]*0.99 > 2*mass_Goldstone:
            # phase_shift = errcl.measurement("phase_shift_%s"%(info["info_string"]), measure_func = calc_phase_shift, sampling_args = ["DONT_RESAMPLE",])
            phase_shift = errcl.measurement("phase_shift_const_%s"%(info["info_string"]), measure_func = calc_phase_shift, sampling_args = ["DONT_RESAMPLE",])
            phase_shift.measure(orig_sample=E_0s, args=(N_L,mass_Goldstone))
            phase_shift.print_to_HDF()



    # ensemble_list = []
    # ensemble_str = "beta=%1.3f, m1/2=%1.3f"
    # color_ind = -1



    # for filename in filelist:
    #     print(filename)
    #     info = HDF_log.get_info_from_HDF5_logfile(filename)
    #     N_T = info[1]
    #     N_L = info[0]
    #     info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(str(info[6]),info[2],info[3],info[4],info[5])
    #     ensemble_str = "beta=%1.3f, m1/2=%1.3f"%(info[3],info[4])
    #     inf_vol = errcl.measurement("infinite_volume_%s_pi"%(info_str), measure_func = None, sampling_args = None)
    #     inf_vol.read_from_HDF()
    #     mass_Goldstone = inf_vol.results["m_inf"].median[0]
    #     energy_levels = errcl.measurement("energy_levels_%s_pipi"%(info[7]), measure_func = None, sampling_args = None)
    #     energy_levels.read_from_HDF()
    #     E_0s = np.swapaxes(energy_levels.results["E_0"].sample,0,1)
    #     E_0s = np.swapaxes(energy_levels.results["E_0"].sample,0,1)
    #     if energy_levels.results["E_0"].median[0]*0.99 > 2*mass_Goldstone:
    #         phase_shift = errcl.measurement("phase_shift_%s"%(info[7]), measure_func = calc_phase_shift, sampling_args = None)
    #         phase_shift.read_from_HDF()
    #     if not (ensemble_str in ensemble_list):
    #         color_ind += 1
    #         ensemble_list.append(ensemble_str)
    #         plt.errorbar(x=phase_shift.results["P_prime"].median, xerr = (phase_shift.results["P_prime"].ep,phase_shift.results["P_prime"].em), y=phase_shift.results["tan_PS"].median, yerr=(phase_shift.results["tan_PS"].ep,phase_shift.results["tan_PS"].em), label = ensemble_str, ls = "", capsize=5, markersize=10, color = color_arr[color_ind])
    #     else:
    #         plt.errorbar(x=phase_shift.results["P_prime"].median, xerr = (phase_shift.results["P_prime"].ep,phase_shift.results["P_prime"].em), y=phase_shift.results["tan_PS"].median, yerr=(phase_shift.results["tan_PS"].ep,phase_shift.results["tan_PS"].em), ls = "", capsize=5, markersize=10, color = color_arr[color_ind])

    # plt.xscale("log")
    # plt.ylabel("tan($\delta$)")
    # plt.xlabel("$\\frac{P}{m_{\pi}}$")
    # plt.legend(fontsize = "xx-small")
    # plt.grid()
    # plt.title(info[7])
    # plt.savefig("plots/tan_PS_"+info[7]+".pdf")
    # plt.show()
    # plt.clf()



    # ensemble_list = []
    # ensemble_str = "beta=%1.3f, m1/2=%1.3f"
    # color_ind = -1
    # for filename in filelist:
    #     print(filename)
    #     info = HDF_log.get_info_from_HDF5_logfile(filename)
    #     N_T = info[1]
    #     N_L = info[0]
    #     info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(str(info[6]),info[2],info[3],info[4],info[5])
    #     ensemble_str = "beta=%1.3f, m1/2=%1.3f"%(info[3],info[4])
    #     inf_vol = errcl.measurement("infinite_volume_%s_pi"%(info_str), measure_func = None, sampling_args = None)
    #     inf_vol.read_from_HDF()
    #     mass_Goldstone = inf_vol.results["m_inf"].median[0]
    #     energy_levels = errcl.measurement("energy_levels_%s_pipi"%(info[7]), measure_func = None, sampling_args = None)
    #     energy_levels.read_from_HDF()
    #     E_0s = np.swapaxes(energy_levels.results["E_0"].sample,0,1)
    #     E_0s = np.swapaxes(energy_levels.results["E_0"].sample,0,1)
    #     if energy_levels.results["E_0"].median[0]*0.99 > 2*mass_Goldstone:
    #         phase_shift = errcl.measurement("phase_shift_%s"%(info[7]), measure_func = calc_phase_shift, sampling_args = None)
    #         phase_shift.read_from_HDF()
    #     if not (ensemble_str in ensemble_list):
    #         color_ind += 1
    #         ensemble_list.append(ensemble_str)
    #         plt.errorbar(x=phase_shift.results["P_2_prime"].median, xerr = (phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em), y=phase_shift.results["P_cot_PS_prime"].median, yerr=(phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em), label = ensemble_str, ls = "", capsize=5, markersize=10, color = color_arr[color_ind])
    #     else:
    #         plt.errorbar(x=phase_shift.results["P_2_prime"].median, xerr = (phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em), y=phase_shift.results["P_cot_PS_prime"].median, yerr=(phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em), ls = "", capsize=5, markersize=10, color = color_arr[color_ind])


    # plt.xscale("log")
    # plt.ylabel("$\\frac{Pcot(\delta)}{m_{\pi}}$")
    # plt.xlabel("$\\frac{P^2}{m_{\pi}}$")
    # plt.legend(fontsize = "xx-small")
    # plt.grid()
    # plt.title(info[7])
    # plt.savefig("plots/P_cot_PS_"+info[7]+".pdf")
    # plt.show()
    # plt.clf()



    # ensemble_list = []
    # ensemble_str = "beta=%1.3f, m1/2=%1.3f"
    # color_ind = -1
    # for filename in filelist:
    #     print(filename)
    #     info = HDF_log.get_info_from_HDF5_logfile(filename)
    #     N_T = info[1]
    #     N_L = info[0]
    #     info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(str(info[6]),info[2],info[3],info[4],info[5])
    #     ensemble_str = "beta=%1.3f, m1/2=%1.3f"%(info[3],info[4])
    #     inf_vol = errcl.measurement("infinite_volume_%s_pi"%(info_str), measure_func = None, sampling_args = None)
    #     inf_vol.read_from_HDF()
    #     mass_Goldstone = inf_vol.results["m_inf"].median[0]
    #     energy_levels = errcl.measurement("energy_levels_%s_pipi"%(info[7]), measure_func = None, sampling_args = None)
    #     energy_levels.read_from_HDF()
    #     E_0s = np.swapaxes(energy_levels.results["E_0"].sample,0,1)
    #     E_0s = np.swapaxes(energy_levels.results["E_0"].sample,0,1)
    #     if energy_levels.results["E_0"].median[0]*0.99 > 2*mass_Goldstone:
    #         phase_shift = errcl.measurement("phase_shift_%s"%(info[7]), measure_func = calc_phase_shift, sampling_args = None)
    #         phase_shift.read_from_HDF()
    #     if not (ensemble_str in ensemble_list):
    #         color_ind += 1
    #         ensemble_list.append(ensemble_str)
    #         plt.errorbar(x=phase_shift.results["P_2_prime"].median, xerr = (phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em), y=phase_shift.results["PS"].median, yerr=(phase_shift.results["PS"].ep,phase_shift.results["PS"].em), label = ensemble_str, ls = "", capsize=5, markersize=10, color = color_arr[color_ind])
    #     else:
    #         plt.errorbar(x=phase_shift.results["P_2_prime"].median, xerr = (phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em), y=phase_shift.results["PS"].median, yerr=(phase_shift.results["PS"].ep,phase_shift.results["PS"].em), ls = "", capsize=5, markersize=10, color = color_arr[color_ind])


    # # plt.xscale("log")
    # plt.ylabel("$\delta$")
    # plt.xlabel("$\\frac{P^2}{m_{\pi}}$")
    # plt.legend(fontsize = "xx-small")
    # plt.grid()
    # plt.title(info[7])
    # plt.savefig("plots/PS_"+info[7]+".pdf")
    # plt.show()
    # plt.clf()

    ############################ MY DATA END ################
        



main()





