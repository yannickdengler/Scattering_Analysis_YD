import numpy as np
import error_classes as errcl
import read_HDF5_logfile as HDF_log 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import generalizedzeta as gz
import plot

# import h5py
# import os

def fit_func_p2(P2, a, b):
    return a + b*P2
def fit_func_p4(P2, a, b, c):
    return a + b*P2 + c*P2*P2
def fit_func_p6(P2, a, b, c, d):
    return a + b*P2 + c*P2*P2 + d*P2*P2*P2

def fit_phase_shift(data, args):       
    result = {}
    print(len(data))
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
    # print("P_2_prime: ", len(P_2_prime), P_2_prime)
    # print("P_cot_PS_prime: ", P_cot_PS_prime)
    result["P_2_prime"] = P_2_prime
    result["P_cot_PS_prime"] = P_cot_PS_prime
    return result

################################ CALCULATION ####################################

def main():
    PATH = "output/result_files/"
    resultfile_list = plot.get_result_files("phase_shift_Fabian")
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
                    # print(len(np.swapaxes(phase_shift.results["P_cot_PS_prime"].sample,0,1)[0]))
                    P_2.append(np.swapaxes(phase_shift.results["P_2_prime"].sample,0,1)[0])
    data = np.zeros((2*len(filenames), len(P_cot_PS[0])))
    for i in range(len(P_cot_PS[0])):
        for j in range(len(filenames)):
            data[j][i] = P_cot_PS[j][i]
    for i in range(len(P_2[0])):
        for j in range(len(filenames)):
            data[j+len(filenames)][i] = P_2[j][i]
    # for i in range(len(filenames)):
    #     print(filenames[i], data[0][i], data[0][i+len(filenames)])
    infos = {}
    infos["filenames"] = filenames
    phase_shift_fit = errcl.measurement("phase_shift_fit", measure_func=fit_phase_shift,sampling_args=("DONT_RESAMPLE",), infos=infos)
    # print(len(data), len(data[0]))
    phase_shift_fit.measure(orig_sample=data, args = None)
    phase_shift_fit.print_to_HDF()
    xarr = np.linspace(0, 1)
    yarr2 = fit_func_p2(xarr, phase_shift_fit.results["a4"].median[0], phase_shift_fit.results["b4"].median[0])
    yarr4 = fit_func_p4(xarr, phase_shift_fit.results["a4"].median[0], phase_shift_fit.results["b4"].median[0], phase_shift_fit.results["c4"].median[0])
    yarr6 = fit_func_p6(xarr, phase_shift_fit.results["a6"].median[0], phase_shift_fit.results["b6"].median[0], phase_shift_fit.results["c6"].median[0], phase_shift_fit.results["d6"].median[0])
    plt.plot(xarr, yarr2)
    plt.plot(xarr, yarr4)
    plt.plot(xarr, yarr6)
    # print(phase_shift_fit.results["P_2"].median)
    plt.scatter(phase_shift_fit.results["P_2_prime"].median, phase_shift_fit.results["P_cot_PS_prime"].median)
    plt.show()


if __name__ == "__main__":
    main()




