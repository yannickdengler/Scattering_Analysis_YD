import os
import error_classes as errcl 
import numpy as np
import matplotlib.pyplot as plt

def write_chi_pt_plot_file():
    m_pi_arr = []
    m_rho_arr = []
    m_pi_rho_arr = []
    f_pi_arr = []
    with open("output/plot_files/chipt_only_stat.dat", "w") as chiptfile:
        for file in os.listdir("output/result_files/"):
            search_str = "phase_shift_fit_P_cot_PS_SP(4)"
            if search_str in file:
                print(file)
                PS_fit_file = errcl.measurement(file[:len(file)-5])
                PS_fit_file.read_from_HDF()
                # a_0 = PS_fit_file.results["a_0_2"].median[0]
                # a_0_p = PS_fit_file.results["a_0_2"].ep[0]
                # a_0_m = PS_fit_file.results["a_0_2"].em[0]
                a_0 = PS_fit_file.results["a_0_0"].median[0]
                a_0_p = PS_fit_file.results["a_0_0"].ep[0]
                a_0_m = PS_fit_file.results["a_0_0"].em[0]
                beta = float(file[35:40])
                m = float(file[43:49])
                m_inf_file = errcl.measurement("infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi"%(beta,m,m))
                m_inf_file.read_from_HDF()
                m_pi = m_inf_file.results["m_inf"].median[0]
                m_pi_e = m_inf_file.results["m_inf"].e[0]
                m_inf_file = errcl.measurement("infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_rho"%(beta,m,m))
                m_inf_file.read_from_HDF()
                m_rho = m_inf_file.results["m_inf"].median[0]
                m_rho_e = m_inf_file.results["m_inf"].e[0]
                # print(a_0, a_0_p, a_0_m, m_pi, m_pi_p, m_pi_m)
                f_pi_data = np.genfromtxt("input/fpi_short.dat")

                for i in range(len(f_pi_data)):
                    if beta == f_pi_data[i][4]:
                        if m == f_pi_data[i][3]:
                            fpi = f_pi_data[i][7]
                            # fpi_e = 1.2*f_pi_data[i][8]
                            fpi_e = f_pi_data[i][8]
                m_f_pi = m_pi/fpi
                m_f_pi_e = np.sqrt((m_pi_e/fpi)**2+(m_pi*fpi_e/fpi**2)**2)
                m_f_pi_2 = (m_pi/fpi)**2
                m_f_pi_2_e = np.sqrt((2*m_pi*m_pi_e/fpi**2)**2+(2*m_pi**2*fpi_e/fpi**3)**2)
                chiptfile.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(a_0, a_0_p, a_0_m, m_f_pi, m_f_pi_e, m_f_pi_2, m_f_pi_2_e, beta, m))
                m_pi_arr.append(m_pi)
                f_pi_arr.append(m_pi/fpi)
                m_pi_rho_arr.append(m_pi/m_rho)

# def write_a_and_b():
#     with open("output/plot_files/a_and_b.dat", "w") as abfile:
#         for file in os.listdir("output/result_files/"):
#             search_str = "phase_shift_fit_P_cot_PS_SP(4)"
#             if search_str in file:
#                 PS_fit_file = errcl.measurement(file[:len(file)-5])
#                 PS_fit_file.read_from_HDF()
#                 a = PS_fit_file.results["a2"].median[0]
#                 a_p = PS_fit_file.results["a2"].ep[0]
#                 a_m = PS_fit_file.results["a2"].em[0]
#                 b = PS_fit_file.results["b2"].median[0]
#                 b_p = PS_fit_file.results["b2"].ep[0]
#                 b_m = PS_fit_file.results["b2"].em[0]
#                 beta = float(file[35:40])
#                 m = float(file[43:49])
#                 m_inf_file = errcl.measurement("infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi"%(beta,m,m))
#                 m_inf_file.read_from_HDF()
#                 mass = m_inf_file.results["m_inf"].median[0]
#                 mass_p = m_inf_file.results["m_inf"].ep[0]
#                 mass_m = m_inf_file.results["m_inf"].em[0]
#                 abfile.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%s\t%e\n"%(a, a_p, a_m, b, b_p, b_m, mass, mass_p, mass_m, beta, m))
    
def write_infinite_volume_extrapolation_data():
    for file in os.listdir("output/result_files/"):
        search_str = "infinite_volume_Fabian"
        if search_str in file:
            inf_vol = errcl.measurement(file[:len(file)-5])
            inf_vol.read_from_HDF()
            gauge_group = inf_vol.infos["gauge_group"]
            if inf_vol.infos["op"] == "pi":
                if inf_vol.infos["level"] == 0:
                    # print(file)
                    # inf_vol.print_everything()
                    N_L = inf_vol.results["N_Ls"].median
                    N_L_inv = inf_vol.results["N_Ls_inv"].median

                    E_pi = inf_vol.results["E"].median
                    E_pi_p = inf_vol.results["E"].ep
                    E_pi_m = inf_vol.results["E"].em
                    m_inf_pi = inf_vol.results["m_inf"].median[0]
                    m_inf_pi_p = inf_vol.results["m_inf"].ep[0]
                    m_inf_pi_m = inf_vol.results["m_inf"].em[0]
                    A_M_pi = inf_vol.results["A_M"].median[0]
                    A_M_pi_p = inf_vol.results["A_M"].ep[0]
                    A_M_pi_m = inf_vol.results["A_M"].em[0]

                    file_rho = file[:len(file)-7]+"rho"
                    inf_vol = errcl.measurement(file_rho)
                    inf_vol.read_from_HDF()
                    E_rho = inf_vol.results["E"].median
                    E_rho_p = inf_vol.results["E"].ep
                    E_rho_m = inf_vol.results["E"].em
                    m_inf_rho = inf_vol.results["m_inf"].median[0]
                    m_inf_rho_p = inf_vol.results["m_inf"].ep[0]
                    m_inf_rho_m = inf_vol.results["m_inf"].em[0]
                    A_M_rho = inf_vol.results["A_M"].median[0]
                    A_M_rho_p = inf_vol.results["A_M"].ep[0]
                    A_M_rho_m = inf_vol.results["A_M"].em[0]

                    file_pipi = file[:len(file)-7]+"pipi"
                    inf_vol = errcl.measurement(file_pipi)
                    inf_vol.read_from_HDF()
                    E_pipi = inf_vol.results["E"].median
                    E_pipi_p = inf_vol.results["E"].ep
                    E_pipi_m = inf_vol.results["E"].em
                    m_inf_pipi = inf_vol.results["m_inf"].median[0]
                    m_inf_pipi_p = inf_vol.results["m_inf"].ep[0]
                    m_inf_pipi_m = inf_vol.results["m_inf"].em[0]
                    A_M_pipi = inf_vol.results["A_M"].median[0]
                    A_M_pipi_p = inf_vol.results["A_M"].ep[0]
                    A_M_pipi_m = inf_vol.results["A_M"].em[0]

                    print(N_L, N_L_inv, E_pi, E_rho, E_pipi, A_M_pi)

                    with open("output/plot_files/infinite_volume/inf_vol_%s_b_%1.3f_m1_%1.3f_m2_%1.3f"%(gauge_group,inf_vol.infos["beta"],inf_vol.infos["m_1"],inf_vol.infos["m_2"]), "w") as filestream:
                        for i in range(len(E_pi)):
                            filestream.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t"%(N_L[i], N_L_inv[i], E_pi[i], E_pi_p[i], E_pi_m[i], m_inf_pi, m_inf_pi_p, m_inf_pi_m, A_M_pi, A_M_pi_p, A_M_pi_m, E_rho[i], E_rho_p[i], E_rho_m[i], m_inf_rho, m_inf_rho_p, m_inf_rho_m, A_M_rho, A_M_rho_p, A_M_rho_m, E_pipi[i], E_pipi_p[i], E_pipi_m[i], m_inf_pipi, m_inf_pipi_p, m_inf_pipi_m, A_M_pipi, A_M_pipi_p, A_M_pipi_m))
                            filestream.write("\n")

def write_phase_shift_fit_data():
    phase_shift_fit_data = []
    ind = 0
    for file in os.listdir("output/result_files/"):
        search_str = "phase_shift_fit_P_cot_PS"
        if search_str in file:
            phase_shift_fit_data.append([])
            phase_shift_fit = errcl.measurement(file[:len(file)-5])
            phase_shift_fit.read_from_HDF()
            phase_shift_fit.print_everything()
            phase_shift_fit_data[ind].append(phase_shift_fit.infos["gauge_group"])
            phase_shift_fit_data[ind].append(phase_shift_fit.infos["beta"])
            phase_shift_fit_data[ind].append(phase_shift_fit.infos["m_1"])
            phase_shift_fit_data[ind].append(phase_shift_fit.infos["m_2"])
            if "a2" in phase_shift_fit.result_names:
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a0"].median[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a0"].ep[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a0"].em[0])
            else:
                phase_shift_fit_data[ind].append(0)
            if "a_0_0" in phase_shift_fit.result_names:
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a_0_0"].median[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a_0_0"].ep[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a_0_0"].em[0])
            else:
                phase_shift_fit_data[ind].append(0)
            if "a2" in phase_shift_fit.result_names:
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a2"].median[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a2"].ep[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a2"].em[0])
            else:
                phase_shift_fit_data[ind].append(0)
            if "b2" in phase_shift_fit.result_names:
                phase_shift_fit_data[ind].append(phase_shift_fit.results["b2"].median[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["b2"].ep[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["b2"].em[0])
            else:
                phase_shift_fit_data[ind].append(0)
            if "a_0_2" in phase_shift_fit.result_names:
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a_0_2"].median[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a_0_2"].ep[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a_0_2"].em[0])
            else:
                phase_shift_fit_data[ind].append(0)
            if "a4" in phase_shift_fit.result_names:
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a4"].median[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a4"].ep[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a4"].em[0])
            else:
                phase_shift_fit_data[ind].append(0)
            if "b4" in phase_shift_fit.result_names:
                phase_shift_fit_data[ind].append(phase_shift_fit.results["b4"].median[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["b4"].ep[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["b4"].em[0])
            else:
                phase_shift_fit_data[ind].append(0)
            if "c4" in phase_shift_fit.result_names:
                phase_shift_fit_data[ind].append(phase_shift_fit.results["c4"].median[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["c4"].ep[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["c4"].em[0])
            else:
                phase_shift_fit_data[ind].append(0)
            if "a_0_4" in phase_shift_fit.result_names:
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a_0_4"].median[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a_0_4"].ep[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a_0_4"].em[0])
            else:
                phase_shift_fit_data[ind].append(0)
            if "a6" in phase_shift_fit.result_names:
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a6"].median[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a6"].ep[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a6"].em[0])
            else:
                phase_shift_fit_data[ind].append(0)
            if "b6" in phase_shift_fit.result_names:
                phase_shift_fit_data[ind].append(phase_shift_fit.results["b6"].median[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["b6"].ep[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["b6"].em[0])
            else:
                phase_shift_fit_data[ind].append(0)
            if "c6" in phase_shift_fit.result_names:
                phase_shift_fit_data[ind].append(phase_shift_fit.results["c6"].median[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["c6"].ep[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["c6"].em[0])
            else:
                phase_shift_fit_data[ind].append(0)
            if "d6" in phase_shift_fit.result_names:
                phase_shift_fit_data[ind].append(phase_shift_fit.results["d6"].median[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["d6"].ep[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["d6"].em[0])
            else:
                phase_shift_fit_data[ind].append(0)
            if "a_0_6" in phase_shift_fit.result_names:
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a_0_6"].median[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a_0_6"].ep[0])
                phase_shift_fit_data[ind].append(phase_shift_fit.results["a_0_6"].em[0])
            else:
                phase_shift_fit_data[ind].append(0)
            ind += 1

        with open("output/plot_files/phase_shift_fit/phase_shift_fit", "w") as filestream:
            for i in range(len(phase_shift_fit_data)):
                filestream.write("%s\t"%phase_shift_fit_data[i][0])
                for j in range(1,len(phase_shift_fit_data[i])):
                    filestream.write("%e\t"%phase_shift_fit_data[i][j])
                filestream.write("\n")

# def write_phase_shift_data():
#     ind = 0
#     for file in os.listdir("output/result_files/"):
#         search_str = "phase_shift_fit_P_cot_PS"
#         if search_str in file:
#             phase_shift_data = []
#             phase_shift_fit = errcl.measurement(file[:len(file)-5])
#             phase_shift_fit.read_from_HDF()







                
                

def main():
    write_chi_pt_plot_file()
    # write_a_and_b()
    # write_infinite_volume_extrapolation_data()
    # write_phase_shift_fit_data()

if __name__ == "__main__":
    main()