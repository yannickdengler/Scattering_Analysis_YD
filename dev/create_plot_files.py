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
                a_0 = PS_fit_file.results["a_0_2"].median[0]
                a_0_p = PS_fit_file.results["a_0_2"].ep[0]
                a_0_m = PS_fit_file.results["a_0_2"].em[0]
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

def write_a_and_b():
    with open("output/plot_files/a_and_b.dat", "w") as abfile:
        for file in os.listdir("output/result_files/"):
            search_str = "phase_shift_fit_P_cot_PS_SP(4)"
            if search_str in file:
                PS_fit_file = errcl.measurement(file[:len(file)-5])
                PS_fit_file.read_from_HDF()
                a = PS_fit_file.results["a2"].median[0]
                a_p = PS_fit_file.results["a2"].ep[0]
                a_m = PS_fit_file.results["a2"].em[0]
                b = PS_fit_file.results["b2"].median[0]
                b_p = PS_fit_file.results["b2"].ep[0]
                b_m = PS_fit_file.results["b2"].em[0]
                beta = float(file[35:40])
                m = float(file[43:49])
                m_inf_file = errcl.measurement("infinite_volume_Fabian_level_0_Scattering_I2_SP(4)_beta%1.3f_m1%1.3f_m2%1.3f_pi"%(beta,m,m))
                m_inf_file.read_from_HDF()
                mass = m_inf_file.results["m_inf"].median[0]
                mass_p = m_inf_file.results["m_inf"].ep[0]
                mass_m = m_inf_file.results["m_inf"].em[0]
                abfile.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%s\t%e\n"%(a, a_p, a_m, b, b_p, b_m, mass, mass_p, mass_m, beta, m))
    
                
                

def main():
    # write_chi_pt_plot_file()
    write_a_and_b()

if __name__ == "__main__":
    main()