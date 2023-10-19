import numpy as np
import error_classes as errcl
from scipy.optimize import curve_fit 


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
    popt, pcov = curve_fit(fit_func_p2, P_2_prime, P_cot_PS_prime)
    a2, b2 = popt
    popt, pcov = curve_fit(fit_func_p4, P_2_prime, P_cot_PS_prime)
    a4, b4, c4 = popt
    popt, pcov = curve_fit(fit_func_p6, P_2_prime, P_cot_PS_prime)
    a6, b6, c6, d6 = popt
    result["a2"], result["b2"], result["a4"], result["b4"], result["c4"], result["a6"], result["b6"], result["c6"], result["d6"] = [[a2,],[b2,],[a4,],[b4,],[c4,],[a6,],[b6,],[c6,],[d6,]] 
    result["a_0_2"] = [1./a2,]
    result["a_0_4"] = [1./a4,]
    result["a_0_6"] = [1./a6,]
    result["P_2_prime"] = P_2_prime
    result["P_cot_PS_prime"] = P_cot_PS_prime
    return result

def fit_func_tan(P2, a, b):
    return a/(a*b*P2+1)

def fit_func_tan_PS(data, args):       
    result = {}
    tan_PS = data[:len(data)//2]
    P_2_prime = data[len(data)//2:]
    popt, pcov = curve_fit(fit_func_tan, P_2_prime, tan_PS)
    a, b = popt
    result["a"] = [a,]
    result["b"] = [b,]
    result["P_2_prime"] = P_2_prime
    result["tan_PS"] = tan_PS
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

def fit_P_cot_PS():
    PATH = "output/result_files/"
    filelist = np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_fit", "str")
    P_cot_PS = []
    P_2 = []
    filenames = []
    for resultfile in resultfile_list:
        phase_shift = errcl.measurement(resultfile)
        if phase_shift.file_exists():
            phase_shift.read_from_HDF()
            P_cot_PS_tmp = np.swapaxes(phase_shift.results["P_cot_PS_prime"].sample,0,1)[0]
            P_2_tmp = np.swapaxes(phase_shift.results["P_2_prime"].sample,0,1)[0]
            if not(any(np.isnan(P_cot_PS_tmp)) or any(np.isinf(P_cot_PS_tmp)) or any(np.isnan(P_2_tmp)) or any(np.isinf(P_2_tmp))):
                if phase_shift.results["P_2_prime"].median[0] < 0.4:
                    filenames.append(resultfile)
                    phase_shift.read_from_HDF()
                    P_cot_PS.append(np.swapaxes(phase_shift.results["P_cot_PS_prime"].sample,0,1)[0])
                    P_2.append(np.swapaxes(phase_shift.results["P_2_prime"].sample,0,1)[0])
    data = np.zeros((2*len(filenames), len(P_cot_PS[0])))
    for i in range(len(P_cot_PS[0])):
        for j in range(len(filenames)):
            data[j][i] = P_cot_PS[j][i]
    for i in range(len(P_2[0])):
        for j in range(len(filenames)):
            data[j+len(filenames)][i] = P_2[j][i]
    infos = {}
    infos["filenames"] = filenames
    phase_shift_fit = errcl.measurement("phase_shift_fit_P_cot_PS", measure_func=fit_func_p_cot_PS,sampling_args=("DONT_RESAMPLE",0), infos=infos)
    phase_shift_fit.measure(orig_sample=data, args = None)
    phase_shift_fit.print_to_HDF()

if __name__ == "__main__":
    # create_all_filenames()
    fit_P_cot_PS()




