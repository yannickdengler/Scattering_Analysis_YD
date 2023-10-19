import numpy as np
import error_classes as errcl
import read_HDF5_logfile as HDF_log 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import plot

color_arr = ["blue", "green", "red", "purple", "lime", "black"]

def infinite_volume(data, args):       
    result = {}
    E_0s = data
    N_Ls = args[0]
    mass_Goldstone_input = args[1]
    N_Ts = args[2]
    mass_Goldstone = 0
    N_Ls_inv = []
    for N_L in N_Ls:
        N_Ls_inv.append(1./N_L)
    def inf_mass_fit_Goldstone(N_L, m_inf, A):
        return m_inf*(1+A*np.exp(-m_inf*N_L)/(m_inf*N_L)**(3/2.))
    def inf_mass_fit_meson(N_L, m_inf, A):
        return m_inf*(1+A*np.exp(-mass_Goldstone_input*N_L)/(mass_Goldstone_input*N_L)**(3/2.))
    if mass_Goldstone_input == 0:
        popt, pcov = curve_fit(f=inf_mass_fit_Goldstone, xdata=N_Ls, ydata=E_0s)
        mass_Goldstone = popt[0]
    else:
        popt, pcov = curve_fit(f=inf_mass_fit_meson, xdata=N_Ls, ydata=E_0s)
        mass_Goldstone = mass_Goldstone_input
    result["mass_Goldstone_input"] = [mass_Goldstone_input,]
    result["mass_Goldstone"] = [mass_Goldstone,]
    result["E"] = E_0s
    result["N_Ts"] = N_Ts
    result["N_Ls"] = N_Ls
    result["N_Ls_inv"] = N_Ls_inv
    result["m_inf"] = [popt[0],]
    result["A_M"] = [popt[1],]
    return result
    
################################ CALCULATION ####################################

def main():
    ops = ("pi", "rho", "pipi")
    filelist_list = []
    with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/HDF5_filelist_infinite_volume", "r") as f:
        for lines in f.readlines():
            if lines[0] != "#":
                filelist_list.append(lines.split())
    for filelist in filelist_list:
        for op in ops:
            E_0 = []
            E_0_const = []
            N_Ls = []
            N_Ts = []
            for files in filelist:
                info = HDF_log.get_info_from_HDF5_logfile(files)
                info["op"] = op
                meas_energlev = errcl.measurement("energy_levels_%s_%s"%(info["info_string"], op))
                meas_energlev.read_from_HDF()
                E_0_const.append(np.swapaxes(meas_energlev.results["E_0_const"].sample,0,1)[0])
                E_0.append(np.swapaxes(meas_energlev.results["E_0"].sample,0,1)[0])
                N_Ls.append(info["N_L"])
                N_Ts.append(info["N_T"])
                if op == "pi":
                    mass_Goldstone = 0
                else:
                    inf_vol = errcl.measurement("infinite_volume_%s_pi"%(info_str))
                    inf_vol.read_from_HDF()
                    mass_Goldstone = inf_vol.results["m_inf"].median[0]
            info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%("I"+str(info["isospin_channel"]),info["gauge_group"],info["beta"],info["m_1"],info["m_2"])
            inf_vol_const = errcl.measurement("infinite_volume_const_%s_%s"%(info_str, op), measure_func = infinite_volume, sampling_args = ("DONT_RESAMPLE",), infos=info)
            inf_vol_const.measure(orig_sample=E_0_const, args=[N_Ls,mass_Goldstone,N_Ts])
            inf_vol_const.print_to_HDF()
            inf_vol = errcl.measurement("infinite_volume_%s_%s"%(info_str, op), measure_func = infinite_volume, sampling_args = ("DONT_RESAMPLE",), infos=info)
            inf_vol.measure(orig_sample=E_0, args=[N_Ls,mass_Goldstone,N_Ts])
            inf_vol.print_to_HDF()

def main_Fabian():
    ops = ("pi", "rho", "pipi")
    beta_m_str = ""
    beta_m_str_list = []
    filelist = plot.get_result_files("energy_levels_Fabian")
    filelist_list = []
    for file in filelist:
        if file[len(file)-5:] == "_pipi":
            print(file)
            meas_energlev = errcl.measurement(file)                     # without hdf5
            meas_energlev.read_from_HDF()
            info = meas_energlev.infos
            beta_m_str = "beta%1.3f_m1%1.3f_m2%1.3f"%(info["beta"],info["m_1"],info["m_2"])
            if not (beta_m_str in beta_m_str_list):
                beta_m_str_list.append(beta_m_str)
                filelist_list.append([])
            filelist_list[beta_m_str_list.index(beta_m_str)].append(file[:len(file)-5])

    for filelist in filelist_list:
        if len(filelist) > 1:
            for op in ops:
                E_0 = []
                E_1 = []
                N_Ls = []
                N_Ts = []
                for files in filelist:
                    print("files: ", files)
                    meas_energlev = errcl.measurement(files+"_"+op)
                    meas_energlev.read_from_HDF()
                    info = meas_energlev.infos                
                    info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%("I"+str(info["isospin_channel"]),info["gauge_group"],info["beta"],info["m_1"],info["m_2"])
                    info["op"] = op
                    info["info_str"] = info_str
                    E_0.append(np.swapaxes(meas_energlev.results["E"].sample,0,1)[0])
                    E_1.append(np.swapaxes(meas_energlev.results["E"].sample,0,1)[1])
                    N_Ls.append(info["N_L"])
                    N_Ts.append(info["N_T"])
                    if op == "pi":
                        mass_Goldstone = 0
                    else:
                        inf_vol = errcl.measurement("infinite_volume_Fabian_level_0_%s_pi"%(info_str))
                        inf_vol.read_from_HDF()
                        mass_Goldstone = inf_vol.results["m_inf"].median[0]
                inf_vol1 = errcl.measurement("infinite_volume_Fabian_level_0_%s_%s"%(info_str, op), measure_func = infinite_volume, sampling_args = ("DONT_RESAMPLE",), infos=info)
                info["level"] = 0
                inf_vol1.measure(orig_sample=E_0, args=[N_Ls,mass_Goldstone,N_Ts])
                inf_vol1.print_to_HDF()
                inf_vol2 = errcl.measurement("infinite_volume_Fabian_level_1_%s_%s"%(info_str, op), measure_func = infinite_volume, sampling_args = ("DONT_RESAMPLE",), infos=info)
                info["level"] = 1
                inf_vol2.measure(orig_sample=E_1, args=[N_Ls,mass_Goldstone,N_Ts])
                inf_vol2.print_to_HDF()




    # with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/HDF5_filelist_infinite_volume", "r") as f:
    #     for lines in f.readlines():
    #         if lines[0] != "#":
    #             filelist_list.append(lines.split())
    # for filelist in filelist_list:
    #     for file in filelist:
    #         print(file)
    #     print()
    # exit()
    # for filelist in filelist_list:
    #     for op in ops:
    #         E_0 = []
    #         E_1 = []
    #         N_Ls = []
    #         for files in filelist:
    #             info = HDF_log.get_info_from_HDF5_logfile(files)
    #             info["op"] = op
    #             meas_energlev = errcl.measurement("energy_levels_Fabian_%s_%s"%(info["info_string"], op))
    #             meas_energlev.read_from_HDF()
    #             E_0.append(np.swapaxes(meas_energlev.results["E"].sample,0,1)[0])
    #             E_1.append(np.swapaxes(meas_energlev.results["E"].sample,0,1)[1])
    #             N_Ls.append(info["N_L"])
    #             if op == "pi":
    #                 mass_Goldstone = 0
    #             else:
    #                 inf_vol = errcl.measurement("infinite_volume_Fabian_%s_pi"%(info_str))
    #                 inf_vol.read_from_HDF()
    #                 mass_Goldstone = inf_vol.results["m_inf"].median[0]
    #         info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%("I"+str(info["isospin_channel"]),info["gauge_group"],info["beta"],info["m_1"],info["m_2"])
    #         inf_vol1 = errcl.measurement("infinite_volume_Fabian_level_1_%s_%s"%(info_str, op), measure_func = infinite_volume, sampling_args = ("DONT_RESAMPLE",), infos=info)
    #         inf_vol1.measure(orig_sample=E_0, args=[N_Ls,mass_Goldstone])
    #         inf_vol1.print_to_HDF()
    #         inf_vol2 = errcl.measurement("infinite_volume_Fabian_level_2_%s_%s"%(info_str, op), measure_func = infinite_volume, sampling_args = ("DONT_RESAMPLE",), infos=info)
    #         inf_vol2.measure(orig_sample=E_1, args=[N_Ls,mass_Goldstone])
    #         inf_vol2.print_to_HDF()













    # ############################ FABIANS DATA ################
    
    # ops = ("pi", "rho", "pipi")
    # PATH = "/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/energy_levels_fabian/"
    # filelist = os.listdir(PATH)
    # # print(filelist)
    # filelist_list = []

    # beta_m_str = ""
    # beta_m_str_list = []
    # for file in filelist:
    #     with h5py.File(PATH+file,"r") as f:
    #         beta_m_str = "beta%1.3f_m1%1.3f_m2%1.3f"%(f["pi/beta"][()],f["pi/m_1"][()],f["pi/m_2"][()])
    #     if not (beta_m_str in beta_m_str_list):
    #         beta_m_str_list.append(beta_m_str)
    #         filelist_list.append([])
    #     # print(beta_m_str_list.index(beta_m_str))
    #     filelist_list[beta_m_str_list.index(beta_m_str)].append(file)

    #                             ########## LEVEL 1 ############

    # for filelist in filelist_list:
    #     if len(filelist) != 1:
    #         for op in ops:
    #             E_0s = []
    #             E_0_errs = []
    #             N_Ls = []
    #             for filename in filelist:
    #                 with h5py.File(PATH+filename,"r") as f:
    #                     # f.visit(print)
    #                     E_0s.append([f[op+"/E"][()][0],])
    #                     E_0_errs.append(f[op+"/Delta_E"][()][0])
    #                     N_Ls.append(f[op+"/N_L"][()])
    #                     if op == "pi":
    #                         mass_Goldstone = 0
    #                     else:
    #                         inf_vol = errcl.measurement("infinite_volume_Fabian_%s_pi"%(info_str))
    #                         inf_vol.read_from_HDF()
    #                         mass_Goldstone = inf_vol.results["m_inf"].median[0]
    #                     info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(f[op+"/isospin_channel"][()],f[op+"/gauge_group"][()].decode(),f[op+"/beta"][()],f[op+"/m_1"][()],f[op+"/m_2"][()])
                
    #             inf_vol = errcl.measurement("infinite_volume_Fabian_%s_%s"%(info_str, op), measure_func = infinite_volume, sampling_args = ("GAUSSIAN",E_0_errs,300,0))
    #             inf_vol.measure(orig_sample=E_0s, args=[N_Ls,mass_Goldstone])
    #             inf_vol.print_to_HDF()

    #                             ########## LEVEL 2 ############

    # for filelist in filelist_list:
    #     if len(filelist) != 1:
    #         for op in ops:
    #             E_0s = []
    #             E_0_errs = []
    #             N_Ls = []
    #             for filename in filelist:
    #                 with h5py.File(PATH+filename,"r") as f:
    #                     # f.visit(print)
    #                     E_0s.append([f[op+"/E"][()][1],])
    #                     E_0_errs.append(f[op+"/Delta_E"][()][1])
    #                     N_Ls.append(f[op+"/N_L"][()])
    #                     if op == "pi":
    #                         mass_Goldstone = 0
    #                     else:
    #                         inf_vol = errcl.measurement("infinite_volume_Fabian_level2_%s_pi"%(info_str))
    #                         inf_vol.read_from_HDF()
    #                         mass_Goldstone = inf_vol.results["m_inf"].median[0]
    #                     info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(f[op+"/isospin_channel"][()],f[op+"/gauge_group"][()].decode(),f[op+"/beta"][()],f[op+"/m_1"][()],f[op+"/m_2"][()])
                
    #             inf_vol = errcl.measurement("infinite_volume_Fabian_level2_%s_%s"%(info_str, op), measure_func = infinite_volume, sampling_args = ("GAUSSIAN",E_0_errs,300,0))
    #             inf_vol.measure(orig_sample=E_0s, args=[N_Ls,mass_Goldstone])
    #             inf_vol.print_to_HDF()


    # filelist_list = []
    # with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/HDF5_filelist_inifinite_volume", "r") as f:
    #     for lines in f.readlines():
    #         if lines[0] != "#":
    #             filelist_list.append(lines.split())
    # for filelist in filelist_list:
    #     if len(filelist) != 1:
    #         for op in ops:
    # # for filelist in filelist_list:
    #             color_ind = 0
    #             info = HDF_log.get_info_from_HDF5_logfile(filelist[0])
    #             E_0s = {}
    #             E_0_eps = {}
    #             E_0_ems = {}
    #             m_inf = {}
    #             m_inf_eps = {}
    #             m_inf_ems = {}
    #             A_M = {}
    #             info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(str(info[6]),info[2],info[3],info[4],info[5])
    #             for op in ops:
    #                 inf_vol = errcl.measurement("infinite_volume_Fabian_%s_%s"%(info_str, op))
    #                 inf_vol.read_from_HDF()
    #                 E_0s[op] = inf_vol.results["E_0s"].median
    #                 E_0_eps[op] = inf_vol.results["E_0s"].ep
    #                 E_0_ems[op] = inf_vol.results["E_0s"].em
    #                 m_inf[op] = inf_vol.results["m_inf"].median[0]
    #                 m_inf_eps[op] = inf_vol.results["m_inf"].ep[0]
    #                 m_inf_ems[op] = inf_vol.results["m_inf"].em[0]
    #                 A_M[op] = inf_vol.results["A_M"].median[0]
    #                 mass_Goldstone = inf_vol.results["mass_Goldstone"].median[0]
    #                 N_Ls = inf_vol.results["N_Ls"].median
    #                 N_Ls_inv = inf_vol.results["N_Ls_inv"].median
    #                 def inf_mass_fit_Goldstone(N_L, m_inf, A):
    #                     return m_inf*(1+A*np.exp(-m_inf*N_L)/(m_inf*N_L)**(3/2.))
    #                 def inf_mass_fit_meson(N_L, m_inf, A):
    #                     return m_inf*(1+A*np.exp(-mass_Goldstone*N_L)/(mass_Goldstone*N_L)**(3/2.))
    #                 xarr_inv = np.linspace(1e-10,max(N_Ls_inv), 100)
    #                 xarr = []
    #                 for x in xarr_inv:
    #                     xarr.append(1./x)
    #                 yarr = []
    #                 yarrp = []
    #                 yarrm = []
    #                 for x in xarr:
    #                     if mass_Goldstone == 0:
    #                         yarr.append(inf_mass_fit_Goldstone(x, m_inf[op], A_M[op]))
    #                         yarrp.append(inf_mass_fit_Goldstone(x, m_inf[op]+m_inf_eps[op], A_M[op]))
    #                         yarrm.append(inf_mass_fit_Goldstone(x, m_inf[op]-m_inf_ems[op], A_M[op]))
    #                     else:
    #                         yarr.append(inf_mass_fit_meson(x, m_inf[op], A_M[op]))
    #                         yarrp.append(inf_mass_fit_meson(x, m_inf[op]+m_inf_eps[op], A_M[op]))
    #                         yarrm.append(inf_mass_fit_meson(x, m_inf[op]-m_inf_ems[op], A_M[op]))
    #                 clr = color_arr[color_ind]
    #                 plt.fill_between(x=xarr_inv, y1=yarrp, y2=yarrm, alpha = 0.5, color=clr)
    #                 plt.errorbar(x=N_Ls_inv, y=E_0s[op], yerr=(E_0_eps[op],E_0_ems[op]), label = op, ls = "", capsize=5, markersize=10, color=clr)
    #                 plt.plot(xarr_inv, yarr, c=clr)
    #                 color_ind += 1
    #             plt.ylabel("E")
    #             plt.xlabel("1/$N_L$")
    #             plt.legend()
    #             plt.grid()
    #             plt.title(info_str+" (Fabian)")
    #             plt.savefig("plots/infinite_volume_Fabian_"+info_str+".pdf")
    #             # plt.show()
    #             plt.clf()

    ############################ FABIANS DATA END ################

    ############################ MY DATA ################


    # ops = ("pi", "rho", "pipi")
    # # ops = ("pi",)
    # filelist_list = []
    # with open("/home/dengler_yannick/Documents/Scattering_Analysis_YD/input/HDF5_filelist_inifinite_volume", "r") as f:
    #     for lines in f.readlines():
    #         if lines[0] != "#":
    #             filelist_list.append(lines.split())
    # for filelist in filelist_list:
    #     for op in ops:
    #         E_0s = []
    #         N_Ls = []
    #         for files in filelist:
    #             print(files)
                # info = HDF_log.get_info_from_HDF5_logfile(files)
    #             meas_energlev = errcl.measurement("energy_levels_%s_%s"%(info[7], op), measure_func = None, sampling_args = (None))
    #             meas_energlev.read_from_HDF()
    #             # E_0s.append(np.swapaxes(meas_energlev.results["E_0"].sample,0,1)[0])
    #             E_0s.append(np.swapaxes(meas_energlev.results["E_0_const"].sample,0,1)[0])
    #             N_Ls.append(info[0])
    #             if op == "pi":
    #                 mass_Goldstone = 0
    #             else:
    #                 inf_vol = errcl.measurement("infinite_volume_%s_pi"%(info_str), measure_func = None, sampling_args = None)
    #                 inf_vol.read_from_HDF()
    #                 mass_Goldstone = inf_vol.results["m_inf"].median[0]
    #         info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(str(info[6]),info[2],info[3],info[4],info[5])
    #         inf_vol = errcl.measurement("infinite_volume_%s_%s"%(info_str, op), measure_func = infinite_volume, sampling_args = ("DONT_RESAMPLE",))
    #         inf_vol.measure(orig_sample=E_0s, args=[N_Ls,mass_Goldstone])
    #         inf_vol.print_to_HDF()

    ############################ MY DATA END ################

    # for filelist in filelist_list:
    #     color_ind = 0
    #     info = HDF_log.get_info_from_HDF5_logfile(filelist[0])
    #     E_0s = {}
    #     E_0_eps = {}
    #     E_0_ems = {}
    #     m_inf = {}
    #     m_inf_eps = {}
    #     m_inf_ems = {}
    #     A_M = {}
    #     info_str = "Scattering_%s_%s_beta%1.3f_m1%1.3f_m2%1.3f"%(str(info[6]),info[2],info[3],info[4],info[5])
    #     for op in ops:
    #         inf_vol = errcl.measurement("infinite_volume_%s_%s"%(info_str, op), measure_func = None, sampling_args = None)
    #         inf_vol.read_from_HDF()
    #         E_0s[op] = inf_vol.results["E_0s"].median
    #         E_0_eps[op] = inf_vol.results["E_0s"].ep
    #         E_0_ems[op] = inf_vol.results["E_0s"].em
    #         m_inf[op] = inf_vol.results["m_inf"].median[0]
    #         m_inf_eps[op] = inf_vol.results["m_inf"].ep[0]
    #         m_inf_ems[op] = inf_vol.results["m_inf"].em[0]
    #         A_M[op] = inf_vol.results["A_M"].median[0]
    #         mass_Goldstone = inf_vol.results["mass_Goldstone"].median[0]
    #         N_Ls = inf_vol.results["N_Ls"].median
    #         N_Ls_inv = inf_vol.results["N_Ls_inv"].median
    #         def inf_mass_fit_Goldstone(N_L, m_inf, A):
    #             return m_inf*(1+A*np.exp(-m_inf*N_L)/(m_inf*N_L)**(3/2.))
    #         def inf_mass_fit_meson(N_L, m_inf, A):
    #             return m_inf*(1+A*np.exp(-mass_Goldstone*N_L)/(mass_Goldstone*N_L)**(3/2.))
    #         xarr_inv = np.linspace(1e-10,max(N_Ls_inv), 100)
    #         xarr = []
    #         for x in xarr_inv:
    #             xarr.append(1./x)
    #         yarr = []
    #         yarrp = []
    #         yarrm = []
    #         for x in xarr:
    #             if mass_Goldstone == 0:
    #                 yarr.append(inf_mass_fit_Goldstone(x, m_inf[op], A_M[op]))
    #                 yarrp.append(inf_mass_fit_Goldstone(x, m_inf[op]+m_inf_eps[op], A_M[op]))
    #                 yarrm.append(inf_mass_fit_Goldstone(x, m_inf[op]-m_inf_ems[op], A_M[op]))
    #             else:
    #                 yarr.append(inf_mass_fit_meson(x, m_inf[op], A_M[op]))
    #                 yarrp.append(inf_mass_fit_meson(x, m_inf[op]+m_inf_eps[op], A_M[op]))
    #                 yarrm.append(inf_mass_fit_meson(x, m_inf[op]-m_inf_ems[op], A_M[op]))
    #         clr = color_arr[color_ind]
    #         plt.fill_between(x=xarr_inv, y1=yarrp, y2=yarrm, alpha = 0.5, color=clr)
    #         plt.errorbar(x=N_Ls_inv, y=E_0s[op], yerr=(E_0_eps[op],E_0_ems[op]), label = op, ls = "", capsize=5, markersize=10, color=clr)
    #         plt.plot(xarr_inv, yarr, c=clr)
    #         color_ind += 1
    #     plt.ylabel("E")
    #     plt.xlabel("1/$N_L$")
    #     plt.legend()
    #     plt.grid()
    #     plt.title(info_str)
    #     # plt.show()
    #     plt.savefig("plots/infinite_volume_"+info_str+".pdf")
    #     plt.clf()



if __name__ == "__main__":
    # main()
    main_Fabian()



