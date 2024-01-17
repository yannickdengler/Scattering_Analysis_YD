import numpy as np
import error_classes as errcl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import os


def fit_func_const(P2, a):
    return a + 0*P2
def fit_func_p2(P2, a, b):
    return a + b*P2
def fit_func_p4(P2, a, b, c):
    return a + b*P2 + c*P2*P2
def fit_func_p6(P2, a, b, c, d):
    return a + b*P2 + c*P2*P2 + d*P2*P2*P2

def fit_func_p_cot_PS(data, args):       
    result = {}
    P_cot_PS_prime = data[:len(data)//2]
    P_2_prime = data[len(data)//2:]
    num = len(P_2_prime)
    popt, pcov = curve_fit(fit_func_const, P_2_prime, P_cot_PS_prime)
    a0 = popt[0]
    result["a0"] = [a0,]
    result["a_0_0"] = [1./a0,]
    if num > 1:
        popt, pcov = curve_fit(fit_func_p2, P_2_prime, P_cot_PS_prime)
        a2, b2 = popt
        result["a2"], result["b2"] = [[a2,],[b2,]]
        result["a_0_2"] = [1./a2,]
        if num > 2:
            popt, pcov = curve_fit(fit_func_p4, P_2_prime, P_cot_PS_prime)
            a4, b4, c4 = popt
            result["a4"], result["b4"], result["c4"] = [[a4,],[b4,],[c4,]]
            result["a_0_4"] = [1./a4,]
            if num > 3:
                popt, pcov = curve_fit(fit_func_p6, P_2_prime, P_cot_PS_prime)
                a6, b6, c6, d6 = popt
                result["a6"], result["b6"], result["c6"], result["d6"] = [[a6,],[b6,],[c6,],[d6,]]
                result["a_0_6"] = [1./a6,]
    result["P_2_prime"] = P_2_prime
    result["P_cot_PS_prime"] = P_cot_PS_prime
    return result

################################ CALCULATION ####################################

def create_all_filenames():
    beta_m_str_list = []
    PATH = "output/result_files/"
    temp = "phase_shift_Fabian"
    filelist = os.listdir(PATH)
    resultfile_list = []
    num = len(temp)
    for file in filelist:
        length = len(file)
        if file[:num] == temp:
            print(file)
            beta_m_str = file[47:74]
            if not(beta_m_str in beta_m_str_list):
                beta_m_str_list.append(beta_m_str)
                resultfile_list.append([])
                resultfile_list[beta_m_str_list.index(beta_m_str)].append(file[:length-5])
            else:
                resultfile_list[beta_m_str_list.index(beta_m_str)].append(file[:length-5])
    for listi in resultfile_list:
        print(listi)
    with open("input/filenames_phase_shift_fit_all", "w") as filestream:
        for filelist in resultfile_list:
            for file in filelist:
                filestream.write(file+"\t")
            filestream.write("\n")

def fit_p_cot_PS(filelist, limit = 0, levels = [0,]):
    P_cot_PS = []
    P_2 = []
    filenames = []
    bm_str = filelist[0][41:65]
    for resultfile in filelist:
        if int(resultfile[25]) in levels:
            phase_shift = errcl.measurement(resultfile)
            if phase_shift.file_exists():
                phase_shift.read_from_HDF()
                if phase_shift.results["P_2_prime"].median[0] < limit or limit == 0:
                    infos = phase_shift.infos
                    infos["levels_fit"] = levels
                    infos["limit_fit"] = limit
                    P_cot_PS_tmp = np.swapaxes(phase_shift.results["P_cot_PS_prime"].sample,0,1)[0]                             # check if you want prime or not
                    P_2_tmp = np.swapaxes(phase_shift.results["P_2_prime"].sample,0,1)[0]
                    if not(any(np.isnan(P_cot_PS_tmp)) or any(np.isinf(P_cot_PS_tmp)) or any(np.isnan(P_2_tmp)) or any(np.isinf(P_2_tmp))):
                        filenames.append(str(resultfile))
                        P_cot_PS.append(np.swapaxes(phase_shift.results["P_cot_PS_prime"].sample,0,1)[0])
                        P_2.append(np.swapaxes(phase_shift.results["P_2_prime"].sample,0,1)[0])
                    else:
                        print("Nans of infs in: "+resultfile)
        
    if len(P_cot_PS) > 1:
        data = np.zeros((2*len(filenames), len(P_cot_PS[0])))
        for i in range(len(P_cot_PS[0])):
            for j in range(len(filenames)):
                data[j][i] = P_cot_PS[j][i]
        for i in range(len(P_2[0])):
            for j in range(len(filenames)):
                data[j+len(filenames)][i] = P_2[j][i]
        for i in range(len(filenames)):
            infos["filenames_fit_%i"%i] = filenames[i]
        phase_shift_fit = errcl.measurement("phase_shift_fit_P_cot_PS_"+bm_str+"_lim_%s"%limit, measure_func=fit_func_p_cot_PS,sampling_args=("DONT_RESAMPLE",0), infos=infos)
        phase_shift_fit.measure(orig_sample=data, args = None)
        phase_shift_fit.print_to_HDF()

def main():
    with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_fit", "r") as f:
        RL = f.readlines()
        for filelist in RL:
            # print(filelist)
            # print()
            fit_p_cot_PS(filelist.split())

if __name__ == "__main__":
    # create_all_filenames()
    main()

