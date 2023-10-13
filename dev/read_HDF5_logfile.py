"""
    @Author: Yannick Dengler
    @Date:   2023-Sep-7
    @Last Modified by: Yannick Dengler
    
    Functions that help with reading the HDF5_logfile
 """

import numpy as np
import h5py

def get_corr_from_HDF5_logfile(filename):
    with h5py.File(filename,"r") as file:
        return file["correlators"][()]

def get_ops_from_HDF5_logfile(filename):
    with h5py.File(filename,"r") as file:
        operators = []
        for ops in file["operators"][()]:
            operators.append(ops.decode())
        return operators

def get_info_from_HDF5_logfile(filename):
    info = {}
    with h5py.File(filename,"r") as file:
        isospin_str = "I%i"%file["isospin_channel"][()]
        info_string = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i"%(isospin_str,file["gauge_group"][()].decode(), file["beta"][()], file["m_1"][()], file["m_2"][()], file["N_T"][()], file["N_L"][()])
        info["N_L"] = file["N_L"][()]
        info["N_T"] = file["N_T"][()]
        info["gauge_group"] = file["gauge_group"][()].decode()
        info["beta"] = file["beta"][()]
        info["m_1"] = file["m_1"][()]
        info["m_2"] = file["m_2"][()]
        info["isospin_channel"] = file["isospin_channel"][()]
        info["info_string"] = info_string
    return info

def get_corr_ops_info_from_HDF5_logfile(filename):
    return (get_corr_from_HDF5_logfile(filename),get_ops_from_HDF5_logfile(filename),get_info_from_HDF5_logfile(filename))

def get_pi_rho_pipi_corr_from_HDF5_logfile(filename):
    corrs = []
    for op in ("pi","rho","pipi"):
        corrs.append(get_corr_from_HDF5_logfile(filename)[get_ops_from_HDF5_logfile(filename).index(op)])
    return corrs

def get_fit_limits(filename):
    info = get_info_from_HDF5_logfile(filename)
    return np.genfromtxt("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/fit_limits/fit_limits_%s"%info["info_string"], int)
