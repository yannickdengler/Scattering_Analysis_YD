"""
    @Author: Yannick Dengler
    @Date:   2023-Sep-8
    @Last Modified by: Yannick Dengler
    
    Takes data and parses it to the resampler and the functions to get an estimation of the error
"""
"""
    Note: Unclear on which error to use. Will continue with BS error (It agrees with std_err like in Gattringer). But they do not agree. Check in with Fabian.
"""
import sys
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import resampling as rs

def mean_orig_sample(orig_sample):
    """
    Takes the mean along the 0th axis. This is used to get the mean result
    """
    mean = []
    for Corr_mont in orig_sample:
        mean.append(np.mean(Corr_mont, axis=0))
    return np.asarray(mean)

class measurement:
    """
    This class contains functions that parses a sample to a resampling algorithm and then to a function that does a calculation
    """
    def __init__(self, name, measure_func = None, sampling_args = ("JK", 0, 0), infos = {}):
        self.name = name
        self.measure_func = measure_func
        self.sampling_args = sampling_args
        self.results = {}
        infos["init"] = True
        self.infos = infos
    def measure(self, orig_sample, args, check_for_bad_results = True):
        """
        Takes a sample and parses it to a function to obtain an estimate for the error
        """
        mean_res = self.measure_func(mean_orig_sample(orig_sample), args)     
        for res_key, res in mean_res.items():
            for val in res:
                if np.isnan(val) or np.isinf(val):
                    sys.exit("NaN or inf in Mean: %s"%(res_key))
        result_samples = {}                                                             
        for key in mean_res:
            result_samples[key] = []
        Resamples = rs.resampling(orig_sample, self.sampling_args)
        tot = len(Resamples)                            
        for i, Resample in zip(range(len(Resamples)), Resamples):
            # print(i*100./tot, "%")
            if check_for_bad_results:
                bad_res = True
                while bad_res:
                    res_tmp = self.measure_func(Resample, args)                               
                    for res in res_tmp.values():
                        if not(any(np.isnan(res)) or any(np.isinf(res))):
                            bad_res = False
            else:
                res_tmp = self.measure_func(Resample, args)                                   
            for key in res_tmp:
                result_samples[key].append(res_tmp[key])
        self.result_names = []
        for key in result_samples:
            self.result_names.append(key)
            self.results[key] = result(key, mean_res[key], result_samples[key], self.sampling_args)
    def visit_print(self):
        with h5py.File("../output/HDF5_logfiles"+self.name+".hdf5","r") as f:
            f.visit(print)
    def print_to_HDF(self, hdfpath = "/home/dengler_yannick/Documents/Scattering_Analysis_YD/output/result_files/"):
        """
        Prints a result to an HDF file
        """
        os.makedirs(hdfpath, exist_ok=True)
        with h5py.File(hdfpath+self.name+".hdf5","w") as f:
            f.create_dataset(str(self.name), data = self.name)
            f.create_dataset("result_names", data = self.result_names)
            for i in range(len(self.sampling_args)):
                f.create_dataset("sampling_args%i"%i, data = self.sampling_args[i])
            for result in self.results.values():
                f.create_dataset(str(result.name)+"_sample", data = result.sample)
                f.create_dataset(str(result.name)+"_result", data = result.result)
            info_names = []
            for key, value in self.infos.items():
                f.create_dataset("info_%s"%key, data = value)
                info_names.append(key)
            f.create_dataset("info_names", data = info_names)
    def read_from_HDF(self, hdfpath = "/home/dengler_yannick/Documents/Scattering_Analysis_YD/output/result_files/"):
        """
        Reads a result from an HDF file
        """
        with h5py.File(hdfpath+self.name+".hdf5","r") as f:
            # f.visit(print)
            # print(f.keys())
            sampling_args_in = []
            self.result_names = []
            for i in range(10):
                tmp = "sampling_args%i"%i
                if tmp in f.keys():
                    tmp2 = f[tmp][()]
                    if isinstance(tmp2, bytes):
                        tmp2 = tmp2.decode()
                    sampling_args_in.append(tmp2)
            self.sampling_args = sampling_args_in
            # print("HEY SMAPLING RGS: ", self.sampling_args)
            for name in f["result_names"]:
                self.result_names.append(name.decode())
            for result_name in self.result_names:
                self.results[result_name] = result(result_name, np.asarray(f[result_name+"_result"]), np.asarray(f[result_name+"_sample"]), self.sampling_args)
            info_names = []
            for name in f["info_names"][()]:
                info_names.append(name.decode())
            for key in info_names:
                val = f["info_%s"%key][()]
                if isinstance(val, bytes):
                    self.infos[key] = f["info_%s"%key][()].decode()
                else:
                    self.infos[key] = f["info_%s"%key][()]
    def print_everything(self):
        print(self.name)
        for key, value in self.infos.items():
            print(key, value)
        for result in self.results.values():
            print(result.name, result.median, result.e)

def JK_err(Sample, res):                                                            
    """
    Calculates the error for a delete-1 Jackknife Sample (not checked)
    """
    print(len(Sample), len(Sample[0]), len(res))

    size = len(res)
    std = np.zeros(size)
    err = np.zeros(size)
    tmp = np.swapaxes(Sample,0,1)
    for i in range(size):
        for val in tmp[i]:
            std[i] += (res[i]-val)**2*((size-1)/size)
            # std[i] += (res[i]-val)**2/size
    for i in range(size):
        err = np.sqrt(std)
    return err

def BS_err(Sample):      
    """
    Calculates the error for a Bootstrap Sample (not checked)
    """               
    num_BS = len(Sample)
    num_res = len(Sample[0])
    var = np.zeros(num_res)    
    err = np.zeros(num_res)
    mean = np.mean(Sample, axis = 0)     
     
    tmp = np.swapaxes(Sample, 0, 1)
    for i in range(num_res):
        for val in tmp[i]:
            var[i] += (mean[i]-val)**2/num_BS
    for i in range(num_res):
        err = np.sqrt(var)
    return err

class result:  
    """
    Class that stores a single result and calculates the error. Is contained in "measurement"
    """          
    def __init__(self, name, result, res_sample, sampling_args = ("JK", 0, 0)):
        self.name = name                                          
        self.sampling_args = sampling_args
        self.result = result
        self.sample = res_sample  
        self.mean = np.mean(res_sample, axis=0)                                                     # I CHANGED THE 0 TO ONE, THIS WAS IMPORTANT
        if not sampling_args[0] == "None":                                # [num_results]
            self.median = []
            self.e = []
            self.ep = []
            self.em = []
            for tmp in np.swapaxes(res_sample,0,1):
                percentage_std = 0.682689
                num = len(tmp)                                                  # num_results
                tmp_2 = np.sort(tmp, kind="mergesort")
                median = tmp_2[num//2]
                self.median.append(median)
                low_ind = math.ceil(num*(1-percentage_std)/2)
                high_ind = math.floor(num*(1+percentage_std)/2)
                ep_tmp = abs(tmp_2[low_ind]-median)
                em_tmp = abs(tmp_2[high_ind]-median)
                self.ep.append(ep_tmp)
                self.em.append(em_tmp)
                self.e.append(max(ep_tmp,em_tmp))
            if sampling_args[0] == "JK" or sampling_args[0] == "JK_SAMEDIM":
                self.e_JK = JK_err(res_sample, result)
            elif sampling_args[0] == "BS_DIFFDIM" or sampling_args[0] == "BS_SAMEDIM":
                self.e_BS = BS_err(res_sample)

    def draw_histogram(self, ind=0, num_bins = 20):  
        """
        Function that draws a histogram from a distribution of a result sample
        """         
        data = np.swapaxes(self.sample,0,1)[ind] 
        start = 0
        stop = len(data)
        counter = 0
        for x in data:
            if x < self.median[ind]+self.ep[ind] and x > self.median[ind]-self.em[ind]:
                counter += 1
        print("%1.3\% of results are in the error (%i\%i)", (counter/(100*stop), counter, stop))
        fig = plt.figure(figsize=(6, 6))
        gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                            left=0.1, right=0.9, bottom=0.1, top=0.9,
                            wspace=0.05, hspace=0.05)
        ax = fig.add_subplot(gs[1, 0])
        ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
        bins = np.linspace(min(data),max(data), num_bins)
        x_axis = []
        for i in range(len(data)):
            x_axis.append(start + i)
        ax.plot(x_axis, data)
        ax.fill_between((0,stop), self.median[ind]-self.em[ind], self.median[ind]+self.ep[ind], color = "red", alpha = 0.5)
        ax.axhline(self.median[ind], color = "red")
        plt.title(self.name)
        plt.hist(data, bins = bins, orientation="horizontal")
        plt.show()

    def percent_nan(self):
        """
        Prints the amount of NaNs in your sample
        """          
        nans = 0
        tot = len(self.sample)
        for res in self.sample:
            if any(np.isnan(res)):
                nans += 1
        print(nans*100./tot,"%")
