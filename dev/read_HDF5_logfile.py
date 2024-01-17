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
        file.visit(print)
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

def get_info_from_Fabian_energy_levels(file):
    info = {}
    isospin_str = "I%i"%file["pi/isospin_channel"][()]
    info_string = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i"%(isospin_str,file["pi/gauge_group"][()].decode(), file["pi/beta"][()], file["pi/m_1"][()], file["pi/m_2"][()], file["pi/N_T"][()], file["pi/N_L"][()])
    info["N_L"] = file["pi/N_L"][()]
    info["N_T"] = file["pi/N_T"][()]
    info["gauge_group"] = file["pi/gauge_group"][()].decode()
    info["beta"] = file["pi/beta"][()]
    info["m_1"] = file["pi/m_1"][()]
    info["m_2"] = file["pi/m_2"][()]
    info["N_hits"] = file["pi/N_hits"][()]
    info["N_mont"] = file["pi/N_mont"][()]
    info["isospin_channel"] = file["pi/isospin_channel"][()]
    info["info_string"] = info_string

    info["Nexp"] = file["pi/Nexp"][()]
    info["antisymmetric"] = file["pi/antisymmetric"][()]
    info["binsize"] = file["pi/binsize"][()]
    info["chi2"] = file["pi/chi2"][()]
    info["dof"] = file["pi/dof"][()]
    info["montecarlotimes"] = file["pi/montecarlotimes"][()]
    info["operators"] = file["pi/operators"][()]
    info["plaquette"] = file["pi/plaquette"][()]
    info["tmax"] = file["pi/tmax"][()]
    info["tmin"] = file["pi/tmin"][()]
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
