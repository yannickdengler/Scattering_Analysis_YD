import numpy as np
import error_classes as errcl
import read_HDF5_logfile as HDF_log 
import matplotlib.pyplot as plt
# from scipy.optimize import curve_fit 
import generalizedzeta as gz

# import h5py
# import os

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
    result["Zeta"] = np.asarray((Zeta,))
    result["M"] = np.asarray((M,))
    result["M_2"] = np.asarray((M*M,))
    return result


################################ CALCULATION ####################################

def main():
    # filelist = os.listdir("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/energy_levels_fabian")
    # print(filelist)
    # for filename in filelist:
    #     with h5py.File(filename,"r") as f:
    #         f.visit(print)





        # info = HDF_log.get_info_from_HDF5_logfile(filename)
        # N_T = info[1]
        # N_L = info[0]
        # info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(str(info[6]),info[2],info[3],info[4],info[5])
        # inf_vol = errcl.measurement("infinite_volume_%s_pi"%(info_str), measure_func = None, sampling_args = None)
        # inf_vol.read_from_HDF()
        # mass_Goldstone = inf_vol.results["m_inf"].median[0]
        # energy_levels = errcl.measurement("energy_levels_%s_pipi"%(info[7]), measure_func = None, sampling_args = None)
        # energy_levels.read_from_HDF()
        # E_0s = np.swapaxes(energy_levels.results["E_0"].sample,0,1)
        # if energy_levels.results["E_0"].median[0]*0.99 > 2*mass_Goldstone:
        #     phase_shift = errcl.measurement("phase_shift_%s"%(info[7]), measure_func = calc_phase_shift, sampling_args = ["DONT_RESAMPLE",])
        #     phase_shift.measure(orig_sample=E_0s, args=(N_L,mass_Goldstone))
        #     phase_shift.print_to_HDF()





    filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/HDF5_filelist_phase_shft", "str")
    for filename in filelist:
        print(filename)
        info = HDF_log.get_info_from_HDF5_logfile(filename)
        N_T = info[1]
        N_L = info[0]
        info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(str(info[6]),info[2],info[3],info[4],info[5])
        inf_vol = errcl.measurement("infinite_volume_%s_pi"%(info_str), measure_func = None, sampling_args = None)
        inf_vol.read_from_HDF()
        mass_Goldstone = inf_vol.results["m_inf"].median[0]
        energy_levels = errcl.measurement("energy_levels_%s_pipi"%(info[7]), measure_func = None, sampling_args = None)
        energy_levels.read_from_HDF()
        E_0s = np.swapaxes(energy_levels.results["E_0"].sample,0,1)
        if energy_levels.results["E_0"].median[0]*0.99 > 2*mass_Goldstone:
            phase_shift = errcl.measurement("phase_shift_%s"%(info[7]), measure_func = calc_phase_shift, sampling_args = ["DONT_RESAMPLE",])
            phase_shift.measure(orig_sample=E_0s, args=(N_L,mass_Goldstone))
            phase_shift.print_to_HDF()

    ensemble_list = []
    ensemble_str = "beta=%1.3f, m1/2=%1.3f"
    color_ind = -1



    for filename in filelist:
        print(filename)
        info = HDF_log.get_info_from_HDF5_logfile(filename)
        N_T = info[1]
        N_L = info[0]
        info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(str(info[6]),info[2],info[3],info[4],info[5])
        ensemble_str = "beta=%1.3f, m1/2=%1.3f"%(info[3],info[4])
        inf_vol = errcl.measurement("infinite_volume_%s_pi"%(info_str), measure_func = None, sampling_args = None)
        inf_vol.read_from_HDF()
        mass_Goldstone = inf_vol.results["m_inf"].median[0]
        energy_levels = errcl.measurement("energy_levels_%s_pipi"%(info[7]), measure_func = None, sampling_args = None)
        energy_levels.read_from_HDF()
        E_0s = np.swapaxes(energy_levels.results["E_0"].sample,0,1)
        E_0s = np.swapaxes(energy_levels.results["E_0"].sample,0,1)
        if energy_levels.results["E_0"].median[0]*0.99 > 2*mass_Goldstone:
            phase_shift = errcl.measurement("phase_shift_%s"%(info[7]), measure_func = calc_phase_shift, sampling_args = None)
            phase_shift.read_from_HDF()
        if not (ensemble_str in ensemble_list):
            color_ind += 1
            ensemble_list.append(ensemble_str)
            plt.errorbar(x=phase_shift.results["P_prime"].median, xerr = (phase_shift.results["P_prime"].ep,phase_shift.results["P_prime"].em), y=phase_shift.results["tan_PS"].median, yerr=(phase_shift.results["tan_PS"].ep,phase_shift.results["tan_PS"].em), label = ensemble_str, ls = "", capsize=5, markersize=10, color = color_arr[color_ind])
        else:
            plt.errorbar(x=phase_shift.results["P_prime"].median, xerr = (phase_shift.results["P_prime"].ep,phase_shift.results["P_prime"].em), y=phase_shift.results["tan_PS"].median, yerr=(phase_shift.results["tan_PS"].ep,phase_shift.results["tan_PS"].em), ls = "", capsize=5, markersize=10, color = color_arr[color_ind])

    plt.xscale("log")
    plt.ylabel("tan($\delta$)")
    plt.xlabel("$\\frac{P}{m_{\pi}}$")
    plt.legend(fontsize = "xx-small")
    plt.grid()
    plt.title(info[7])
    plt.savefig("plots/tan_PS_"+info[7]+".pdf")
    plt.show()
    plt.clf()



    ensemble_list = []
    ensemble_str = "beta=%1.3f, m1/2=%1.3f"
    color_ind = -1
    for filename in filelist:
        print(filename)
        info = HDF_log.get_info_from_HDF5_logfile(filename)
        N_T = info[1]
        N_L = info[0]
        info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(str(info[6]),info[2],info[3],info[4],info[5])
        ensemble_str = "beta=%1.3f, m1/2=%1.3f"%(info[3],info[4])
        inf_vol = errcl.measurement("infinite_volume_%s_pi"%(info_str), measure_func = None, sampling_args = None)
        inf_vol.read_from_HDF()
        mass_Goldstone = inf_vol.results["m_inf"].median[0]
        energy_levels = errcl.measurement("energy_levels_%s_pipi"%(info[7]), measure_func = None, sampling_args = None)
        energy_levels.read_from_HDF()
        E_0s = np.swapaxes(energy_levels.results["E_0"].sample,0,1)
        E_0s = np.swapaxes(energy_levels.results["E_0"].sample,0,1)
        if energy_levels.results["E_0"].median[0]*0.99 > 2*mass_Goldstone:
            phase_shift = errcl.measurement("phase_shift_%s"%(info[7]), measure_func = calc_phase_shift, sampling_args = None)
            phase_shift.read_from_HDF()
        if not (ensemble_str in ensemble_list):
            color_ind += 1
            ensemble_list.append(ensemble_str)
            plt.errorbar(x=phase_shift.results["P_2_prime"].median, xerr = (phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em), y=phase_shift.results["P_cot_PS_prime"].median, yerr=(phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em), label = ensemble_str, ls = "", capsize=5, markersize=10, color = color_arr[color_ind])
        else:
            plt.errorbar(x=phase_shift.results["P_2_prime"].median, xerr = (phase_shift.results["P_2_prime"].ep,phase_shift.results["P_2_prime"].em), y=phase_shift.results["P_cot_PS_prime"].median, yerr=(phase_shift.results["P_cot_PS_prime"].ep,phase_shift.results["P_cot_PS_prime"].em), ls = "", capsize=5, markersize=10, color = color_arr[color_ind])


    plt.xscale("log")
    plt.ylabel("$\\frac{Pcot(\delta)}{m_{\pi}}$")
    plt.xlabel("$\\frac{P^2}{m_{\pi}}$")
    plt.legend(fontsize = "xx-small")
    plt.grid()
    plt.title(info[7])
    plt.savefig("plots/P_cot_PS_"+info[7]+".pdf")
    plt.show()
    plt.clf()



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
        



main()





