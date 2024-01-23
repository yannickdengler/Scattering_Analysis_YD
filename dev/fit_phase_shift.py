import numpy as np
import error_classes as errcl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import os
from scipy.integrate import quad as integrate

def get_max_min_int(range_x, coeffs):
    xarr = np.linspace(range_x[0], range_x[1])
    yarr = np.zeros(len(xarr))
    if len(coeffs) == 1:
        integral = integrate(fit_func_const, range_x[0], range_x[1], args = tuple(coeffs))[0]
        for i in range(len(xarr)):
            yarr[i] = fit_func_const(xarr[i], coeffs[0])
    elif len(coeffs) == 2:
        integral = integrate(func=fit_func_p2, a=range_x[0], b=range_x[1], args = tuple(coeffs))[0]
        for i in range(len(xarr)):
            yarr[i] = fit_func_p2(xarr[i], coeffs[0], coeffs[1])
    elif len(coeffs) == 3:
        integral = integrate(fit_func_p4, range_x[0], range_x[1], args = tuple(coeffs))[0]
        for i in range(len(xarr)):
            yarr[i] = fit_func_p4(xarr[i], coeffs[0], coeffs[1], coeffs[2])
    elif len(coeffs) == 4:
        integral = integrate(fit_func_p6, range_x[0], range_x[1], args = tuple(coeffs))[0]
        for i in range(len(xarr)):
            yarr[i] = fit_func_p6(xarr[i], coeffs[0], coeffs[1], coeffs[2], coeffs[3])
    else:
         return None
    minval = min(yarr)
    maxval = max(yarr)
    return integral, minval, maxval, (minval+maxval)/2

def get_indices_from_index(index):
    indices = [0,0]
    for i in range(index):
        if indices[1]+1 > indices[0]:
            indices[0] += 1
            indices[1] = 0
        else:
            indices[1]+=1
    return indices

# for i in range(20):
#     print(i, get_indices_from_index(i))
# exit()

def get_mean_difference(data, args):                    #  data - bootstrap of coeffs, args - array of coeffs [order]x[a,b,c, ...]
    result = {}                                                            
    limit = args[1]
    coeffs_sample = np.zeros((10,10))
    mean_coeffs = np.zeros((10,10))
    for i in range(len(data)):
        indices = get_indices_from_index(i)
        coeffs_sample[indices[0]][indices[1]] = data[i]
        mean_coeffs[indices[0]][indices[1]] = args[0][i]
    coeffs_tmp = []
    indices = get_indices_from_index(len(data))
    for i in range(indices[0]):
        coeffs_tmp.append([])
        for j in range(i+1):
            coeffs_tmp[i].append(coeffs_sample[i][j]-mean_coeffs[i][j])

    for i in range(len(coeffs_tmp)):
        if i == 0:
            result["integral"+str(2*i)] = [integrate(fit_func_const, 0, limit, args = tuple(coeffs_tmp[i]))[0],]
        elif i == 1:
            result["integral"+str(2*i)] = [integrate(fit_func_p2, 0, limit, args = tuple(coeffs_tmp[i]))[0],]
        elif i == 2:
            result["integral"+str(2*i)] = [integrate(fit_func_p4, 0, limit, args = tuple(coeffs_tmp[i]))[0],]
        elif i == 3:
            result["integral"+str(2*i)] = [integrate(fit_func_p6, 0, limit, args = tuple(coeffs_tmp[i]))[0],]
    return result


        


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
    limit = args[0]
    range_x = [0, limit]
    P_cot_PS_prime = data[:len(data)//2]
    P_2_prime = data[len(data)//2:]
    num = len(P_2_prime)
    popt, pcov = curve_fit(fit_func_const, P_2_prime, P_cot_PS_prime)
    a0 = popt[0]
    # integral0, minval0, maxval0, minmaxval0 = get_max_min_int(range_x=[min(P_2_prime), max(P_2_prime)], coeffs=[a0,])
    # integral0, minval0, maxval0, minmaxval0 = get_max_min_int(range_x=range_x, coeffs=[a0,])
    result["a0"] = [a0,]
    # result["integral0"], result["minval0"], result["maxval0"], result["minmaxval0"] = [[integral0,],[minval0,],[maxval0,],[minmaxval0,]]
    result["a_0_0"] = [1./a0,]
    if num > 2:
        popt, pcov = curve_fit(fit_func_p2, P_2_prime, P_cot_PS_prime)
        a2, b2 = popt
        # integral2, minval2, maxval2, minmaxval2 = get_max_min_int(range_x=[min(P_2_prime), max(P_2_prime)], coeffs=[a2,b2])
        # integral2, minval2, maxval2, minmaxval2 = get_max_min_int(range_x=range_x, coeffs=[a0,])
        result["a2"], result["b2"] = [[a2,],[b2,]]
        # result["integral2"], result["minval2"], result["maxval2"], result["minmaxval2"] = [[integral2,],[minval2,],[maxval2,],[minmaxval2,]]
        result["a_0_2"] = [1./a2,]
        if num > 3:
            popt, pcov = curve_fit(fit_func_p4, P_2_prime, P_cot_PS_prime)
            a4, b4, c4 = popt
            # integral4, minval4, maxval4, minmaxval4 = get_max_min_int(range_x=[min(P_2_prime), max(P_2_prime)], coeffs=[a4,b4,c4])
            # integral4, minval4, maxval4, minmaxval4 = get_max_min_int(range_x=range_x, coeffs=[a0,])
            result["a4"], result["b4"], result["c4"] = [[a4,],[b4,],[c4,]]
            # result["integral4"], result["minval4"], result["maxval4"], result["minmaxval4"] = [[integral4,],[minval4,],[maxval4,],[minmaxval4,]]
            result["a_0_4"] = [1./a4,]
            if num > 4:
                popt, pcov = curve_fit(fit_func_p6, P_2_prime, P_cot_PS_prime)
                a6, b6, c6, d6 = popt
                # integral6, minval6, maxval6, minmaxval6 = get_max_min_int(range_x=[min(P_2_prime), max(P_2_prime)], coeffs=[a6,b6,c6,d6])
                # integral6, minval6, maxval6, minmaxval6 = get_max_min_int(range_x=range_x, coeffs=[a0,])
                result["a6"], result["b6"], result["c6"], result["d6"] = [[a6,],[b6,],[c6,],[d6,]]
                # result["integral6"], result["minval6"], result["maxval6"], result["minmaxval6"] = [[integral6,],[minval6,],[maxval6,],[minmaxval6,]]
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
        phase_shift_fit.measure(orig_sample=data, args = [limit,])                
        phase_shift_fit.add_single_result(phase_shift.results["mass_Goldstone"])
        # phase_shift_fit.print_everything()
        # print(len(data))
        # print(len(data[0]))

        # print(phase_shift_fit.result_names)

        result_name_arr = ["a0", "a2", "b2", "a4", "b4", "c4", "a6", "b6", "c6", "d6"]
        integral_data = []
        mean_coeffs = []
        for i in range(len(result_name_arr)):
            if result_name_arr[i] in phase_shift_fit.result_names:
                integral_data.append(np.transpose(phase_shift_fit.results[result_name_arr[i]].sample)[0])
                mean_coeffs.append(phase_shift_fit.results[result_name_arr[i]].median[0])

        integral_meas = errcl.measurement("integral_data", measure_func=get_mean_difference, sampling_args=("DONT_RESAMPLE",0), infos = phase_shift_fit.infos)
        integral_meas.measure(orig_sample=integral_data, args = [mean_coeffs, limit])
        for key in integral_meas.result_names:
            phase_shift_fit.add_single_result(integral_meas.results[key])


        # exit()


        phase_shift_fit.print_to_HDF()
        # phase_shift_fit.print_everything()

def main():
    with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/filenames_phase_shift_fit", "r") as f:
        RL = f.readlines()
        for filelist in RL:
            # print(filelist)
            # print()
            fit_p_cot_PS(filelist.split(), limit = 0.9)

if __name__ == "__main__":
    # create_all_filenames()
    main()

