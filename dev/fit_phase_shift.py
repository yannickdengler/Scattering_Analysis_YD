import numpy as np
import error_classes as errcl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import os


def fit_func_p2(P2, a, b):
    return a + b*P2
def fit_func_p4(P2, a, b, c):
    return a + b*P2 + c*P2*P2
def fit_func_p6(P2, a, b, c, d):
    return a + b*P2 + c*P2*P2 + d*P2*P2*P2
def fit_func_p4_fixed_400(P2, b, c):
    return -1/400 + b*P2 + c*P2*P2
def fit_func_p2_fixed_400(P2, b):
    return -1/400 + b*P2

def fit_func_p_cot_PS(data, args):       
    result = {}
    P_cot_PS_prime = data[:len(data)//2]
    P_2_prime = data[len(data)//2:]
    num = len(P_2_prime)
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
    PATH = "output/result_files/"
    temp = "phase_shift_Fabian"
    filelist = os.listdir(PATH)
    resultfile_list = []
    num = len(temp)
    for file in filelist:
        length = len(file)
        if file[:num] == temp:
            resultfile_list.append(file[:length-5])         

    with open("input/filenames_phase_shift_fit_all", "w") as file:
        for filename in resultfile_list:
            file.write("%s"%filename+"\n")

def fit_p_cot_PS(filelist, limit = 0):
    PATH = "output/result_files/"
    # filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_fit", "str")
    P_cot_PS = []
    P_2 = []
    filenames = []
    bm_str = filelist[0][41:65]
    for resultfile in filelist:
        phase_shift = errcl.measurement(resultfile)
        print(resultfile)
        if phase_shift.file_exists():
            phase_shift.read_from_HDF()
            P_cot_PS_tmp = np.swapaxes(phase_shift.results["P_cot_PS_prime"].sample,0,1)[0]
            P_2_tmp = np.swapaxes(phase_shift.results["P_2_prime"].sample,0,1)[0]
            # P_cot_PS_tmp = np.swapaxes(phase_shift.results["P_cot_PS"].sample,0,1)[0]
            # P_2_tmp = np.swapaxes(phase_shift.results["P_2"].sample,0,1)[0]
            print(not(any(np.isnan(P_cot_PS_tmp)) or any(np.isinf(P_cot_PS_tmp)) or any(np.isnan(P_2_tmp)) or any(np.isinf(P_2_tmp))))
            if not(any(np.isnan(P_cot_PS_tmp)) or any(np.isinf(P_cot_PS_tmp)) or any(np.isnan(P_2_tmp)) or any(np.isinf(P_2_tmp))):
                filenames.append(str(resultfile))
                phase_shift.read_from_HDF()
                P_cot_PS.append(np.swapaxes(phase_shift.results["P_cot_PS_prime"].sample,0,1)[0])
                P_2.append(np.swapaxes(phase_shift.results["P_2_prime"].sample,0,1)[0])
                # P_cot_PS.append(np.swapaxes(phase_shift.results["P_cot_PS"].sample,0,1)[0])
                # P_2.append(np.swapaxes(phase_shift.results["P_2"].sample,0,1)[0])
    num = len(P_cot_PS)
    if num > 1:
        data = np.zeros((2*len(filenames), len(P_cot_PS[0])))
        for i in range(len(P_cot_PS[0])):
            for j in range(len(filenames)):
                data[j][i] = P_cot_PS[j][i]
        for i in range(len(P_2[0])):
            for j in range(len(filenames)):
                data[j+len(filenames)][i] = P_2[j][i]
        infos = {}
        for i in range(len(filenames)):
            infos["filenames_%i"%i] = filenames[i]
        phase_shift_fit = errcl.measurement("phase_shift_fit_P_cot_PS_"+bm_str, measure_func=fit_func_p_cot_PS,sampling_args=("DONT_RESAMPLE",0), infos=infos)
        phase_shift_fit.measure(orig_sample=data, args = None)
        for resultfile in filelist:
            phase_shift = errcl.measurement(resultfile)
            if phase_shift.file_exists():
                phase_shift.read_from_HDF()
                phase_shift_fit.add_all_results(phase_shift)
        phase_shift_fit.print_everything()
        phase_shift_fit.print_to_HDF()



def fit_all_p_cot_PS(filelist, limit = 0):
    P_cot_PS = []
    P_2 = []
    filenames = []
    for lines in filelist:
        words = lines.split()
        bm_str = words[0][41:65]
        for resultfile in words:
            phase_shift = errcl.measurement(resultfile)
            if phase_shift.file_exists():
                phase_shift.read_from_HDF()
                P_cot_PS_tmp = np.swapaxes(phase_shift.results["P_cot_PS_prime"].sample,0,1)[0]
                # P_cot_PS_tmp = np.swapaxes(phase_shift.results["P_cot_PS"].sample,0,1)[0]
                P_2_tmp = np.swapaxes(phase_shift.results["P_2_prime"].sample,0,1)[0]
                # P_2_tmp = np.swapaxes(phase_shift.results["P_2"].sample,0,1)[0]
                # print(not(any(np.isnan(P_cot_PS_tmp)) or any(np.isinf(P_cot_PS_tmp)) or any(np.isnan(P_2_tmp)) or any(np.isinf(P_2_tmp))))
                if not(any(np.isnan(P_cot_PS_tmp)) or any(np.isinf(P_cot_PS_tmp)) or any(np.isnan(P_2_tmp)) or any(np.isinf(P_2_tmp))):
                    # print(phase_shift.results["P_2_prime"].median[0], limit)
                    if phase_shift.results["P_2_prime"].median[0] < limit:
                    # print(phase_shift.results["P_2"].median[0], limit)
                    # if phase_shift.results["P_2"].median[0] < limit:
                        print(phase_shift.name)
                        filenames.append(str(resultfile))
                        P_cot_PS.append(np.swapaxes(phase_shift.results["P_cot_PS_prime"].sample,0,1)[0])
                        P_2.append(np.swapaxes(phase_shift.results["P_2_prime"].sample,0,1)[0])
                        # P_cot_PS.append(np.swapaxes(phase_shift.results["P_cot_PS"].sample,0,1)[0])
                        # P_2.append(np.swapaxes(phase_shift.results["P_2"].sample,0,1)[0])
    num = len(P_cot_PS)
    print(num)
    if num > 1:
        data = np.zeros((2*len(filenames), len(P_cot_PS[0])))
        for i in range(len(P_cot_PS[0])):
            for j in range(len(filenames)):
                data[j][i] = P_cot_PS[j][i]
        for i in range(len(P_2[0])):
            for j in range(len(filenames)):
                data[j+len(filenames)][i] = P_2[j][i]
        infos = {}
        for i in range(len(filenames)):
            infos["filenames_%i"%i] = filenames[i]
        infos["limit"] = limit
        # phase_shift_fit = errcl.measurement("phase_shift_fit_P_cot_PS_lim_%e_"%limit+bm_str, measure_func=fit_func_p_cot_PS,sampling_args=("DONT_RESAMPLE",0), infos=infos)
        phase_shift_fit = errcl.measurement("phase_shift_fit_P_cot_PS_lim_%e_no_prime_"%limit+bm_str, measure_func=fit_func_p_cot_PS,sampling_args=("DONT_RESAMPLE",0), infos=infos)
        phase_shift_fit.measure(orig_sample=data, args = None)
        for resultfile in filelist:
            phase_shift = errcl.measurement(resultfile)
            if phase_shift.file_exists():
                phase_shift.read_from_HDF()
                phase_shift_fit.add_all_results(phase_shift)
        # phase_shift_fit.print_everything()
        phase_shift_fit.print_to_HDF()



if __name__ == "__main__":
    create_all_filenames()
    with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift", "r") as f:
        RL = f.readlines()    
        print(RL)    
        fit_p_cot_PS(RL, limit=1e50)






    # with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_plot_SP4_lev0_split", "r") as f:
    #     RL = f.readlines()
        # for lines in f.readlines():
        #     filelist = lines.split()
        #     fit_p_cot_PS(filelist)
        # for lim in np.linspace(start=2e-2,stop=0.5, num = 5):
        #     print(lim)
        #     fit_all_p_cot_PS(RL, limit=lim)
        # fit_all_p_cot_PS(RL, limit=1e50)


    





    




