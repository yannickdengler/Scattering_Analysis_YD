import numpy as np
import error_classes as errcl
import generalizedzeta as gz
import plot

import os

color_arr = ["blue", "green", "red", "purple", "lime", "black", "grey", "orange", "olive", "skyblue", "fuchsia", "peru", "firebrick"]


def generalized_mom(E_pipi, m_meson):
    if E_pipi < 2*m_meson:
        return 0
    return 2*np.arcsin(np.sqrt(0.5*(np.cosh(E_pipi/2.)-np.cosh(m_meson))))

def generalized_mom_non_lattice(E_pipi, m_meson):
    if E_pipi < 2*m_meson:
        return 0
    return np.sqrt(E_pipi*E_pipi/4-m_meson*m_meson)

def calc_phase_shift(data, args):       
    result = {}
    E_pipi = data[0]
    N_L = args[0]
    mass_Goldstone = args[1]
    N_L_inv = 1./N_L
    # P = generalized_mom_non_lattice(E_pipi, mass_Goldstone)
    P = generalized_mom(E_pipi, mass_Goldstone)
    if P == 0:
        P = 1e-99
        q = 1e-99
        Zeta = -1e99
        tan_PS = 1e-99
    else:
        q = P*N_L/(2*np.pi)
        Zeta = gz.Zeta(q**2)
        tan_PS = np.pi**(3/2.)*q/Zeta

    M = 16*np.pi*tan_PS/np.sqrt(tan_PS**2+1)                                                # sin(arctan(tan(x)))
    result["N_L"] = np.asarray((N_L,))
    result["N_L_inv"] = np.asarray((N_L_inv,))
    result["mass_Goldstone"] = np.asarray((mass_Goldstone,))
    result["E_pipi"] = np.asarray((E_pipi,))
    result["E_pipi_prime"] = np.asarray((E_pipi/mass_Goldstone,))
    result["E_m_2m"] = np.asarray((E_pipi-2*mass_Goldstone,))
    result["E_m_2m_prime"] = np.asarray((E_pipi/mass_Goldstone-2,))
    result["s_pipi"] = np.asarray((E_pipi**2,))
    result["s_pipi_prime"] = np.asarray(((E_pipi/mass_Goldstone)**2,))
    result["P"] = np.asarray((P,))
    result["P_prime"] = np.asarray((P/mass_Goldstone,))
    result["P_2"] = np.asarray((P**2,))
    result["P_2_prime"] = np.asarray(((P/mass_Goldstone)**2,))
    result["q"] = np.asarray((q,))
    result["q_2"] = np.asarray((q**2,))
    result["tan_PS"] = np.asarray((tan_PS,))
    result["tan_PS/P"] = np.asarray((tan_PS/P,))
    result["tan_PS/P_prime"] = np.asarray((tan_PS*mass_Goldstone/P,))
    result["PS"] = np.asarray((np.arctan(tan_PS),))
    result["cot_PS"] = np.asarray((1./tan_PS,))
    result["P_cot_PS"] = np.asarray((P/tan_PS,))
    result["P_cot_PS_prime"] = np.asarray((P/(tan_PS*mass_Goldstone),))
    result["P3_cot_PS/E"] = np.asarray((P**3/(E_pipi*tan_PS*mass_Goldstone**2),))
    result["P3_cot_PS/E_prime"] = np.asarray((P/(E_pipi*tan_PS*mass_Goldstone**2),))
    result["Zeta"] = np.asarray((Zeta,))
    result["M"] = np.asarray((M,))
    result["M_2"] = np.asarray((M*M,))
    # result["sigma"] = M**2/(16*np.pi*E_pipi**2)
    # result["sigma_prime"] = M**2*mass_Goldstone**2/(16*np.pi*E_pipi**2)
    # result["velocity"] = P/np.sqrt(mass_Goldstone**2+P**2)
    if np.isnan(P/(tan_PS*mass_Goldstone)) or np.isinf(P/(tan_PS*mass_Goldstone)):
        print("NaN or inf!")
        print(P, q, Zeta, tan_PS, P/(tan_PS*mass_Goldstone), mass_Goldstone)
    return result

################################ CALCULATION ####################################

def create_all_filenames():
    PATH = "output/result_files/"
    temp = "infinite_volume_Fabian"
    filelist = os.listdir(PATH)
    resultfile_list = []
    num = len(temp)
    for file in filelist:
        length = len(file)
        if file[:num] == temp:
            if file[length-10:] == "_pipi.hdf5":
                resultfile_list.append(file[:length-5])         

    with open("input/filenames_phase_shift_all", "w") as file:
        for filename in resultfile_list:
            file.write("%s"%filename+"\n")

def main():
    PATH = "output/result_files/"
    filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_all", "str")

    for resultfile in resultfile_list:
        if resultfile[len(resultfile)-5:] == "_pipi":
            energy_levels = errcl.measurement(resultfile)
            energy_levels.read_from_HDF()
            info = energy_levels.infos
            info_str_inf = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%("I"+str(info["isospin_channel"]),info["gauge_group"],info["beta"],info["m_1"],info["m_2"])
            if info["gauge_group"] != "SU(3)":
                if ("infinite_volume_Fabian_level_1_%s_pi"%(info_str_inf)+".hdf5") in os.listdir(PATH):
                    inf_vol_rho = errcl.measurement("infinite_volume_Fabian_level_1_%s_rho"%(info_str_inf))
                    inf_vol_rho.read_from_HDF()
                    inf_vol_pi = errcl.measurement("infinite_volume_Fabian_level_1_%s_pi"%(info_str_inf))
                    inf_vol_pi.read_from_HDF()
                    mass_Goldstone = inf_vol_pi.results["m_inf"].median[0]
                    m_rho_m_pi = inf_vol_rho.results["m_inf"].median[0]/inf_vol_pi.results["m_inf"].median[0]
                    for i in range(len(energy_levels.results["E"].median)):
                        if (energy_levels.results["E"].median[i]*0.99 > 2*mass_Goldstone) and (energy_levels.results["E"].median[i]*1.01 < 4*mass_Goldstone):
                            E = [np.swapaxes(energy_levels.results["E"].sample,0,1)[i],]
                            phase_shift = errcl.measurement("phase_shift_Fabian_%i_%s"%(i, info["info_string"]), measure_func = calc_phase_shift, sampling_args = ["DONT_RESAMPLE",50], infos=info)
                            phase_shift.measure(orig_sample=E, args=(info["N_L"],mass_Goldstone))
                            ########################### CHECK N_L
                            phase_shift.infos["m_rho_m_pi"] = m_rho_m_pi
                            phase_shift.print_to_HDF()

def main_Fabian():
    filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_all", "str")
    # filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_SU3", "str")

    for file in filelist:
        # print(file)
        inf_vol = errcl.measurement(file)
        inf_vol.read_from_HDF()
        N_Ls = inf_vol.results["N_Ls"].median
        info = inf_vol.infos
        for i in range(len(inf_vol.results["E"].median)):  
            info_string = "Scattering_I%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i"%(info["isospin_channel"],info["gauge_group"], info["beta"],info["m_1"],info["m_2"], inf_vol.results["N_Ts"].median[i], inf_vol.results["N_Ls"].median[i])
            mass_Goldstone = inf_vol.results["mass_Goldstone"].median[0]
            if (inf_vol.results["E"].median[i]-inf_vol.results["E"].em[i] > 2*mass_Goldstone) and (inf_vol.results["E"].median[i]+inf_vol.results["E"].ep[i] < 4*mass_Goldstone):
                phase_shift = errcl.measurement("phase_shift_Fabian_level_%i_%s"%(info["level"], info_string), measure_func = calc_phase_shift, sampling_args = ["DONT_RESAMPLE",0], infos=info)
                # phase_shift = errcl.measurement("phase_shift_Fabiantest_level_%i_%s"%(info["level"], info_string), measure_func = calc_phase_shift, sampling_args = ["None",0], infos=info)
                phase_shift.measure(orig_sample=[np.swapaxes(inf_vol.results["E"].sample,0,1)[i],], args=(N_Ls[i],mass_Goldstone))
                phase_shift.add_single_result(inf_vol.results["mass_Goldstone"])
                phase_shift.print_to_HDF()
                # phase_shift.print_everything()
                # exit()

def print_phase_shift():
    filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_all", "str")
    print(filelist)
    for file in filelist:
        # print(file)
        inf_vol = errcl.measurement(file)
        inf_vol.read_from_HDF()
        print(file,inf_vol.infos["beta"],inf_vol.infos["m_1"],inf_vol.infos["N_mont"])
        # print(file, "\t\t\t", inf_vol.infos["N_mont"])
        # inf_vol.print_everything()



if __name__ == "__main__":
    # create_all_filenames()
    # main()
    # main_Fabian()
    print_phase_shift()




