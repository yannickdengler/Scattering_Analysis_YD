def calc_corr(Correlators, args):                               # Correlators [Ops][N_T] 
    Op_names = args[0]                                                       # Initialization of args   # same length as Correlators
    if len(Correlators) != len(Op_names):
        print(len(Correlators),len(Op_names))
        print("length of Corr is not same as Ops, calc_corr")
        return 1
    result = {}                                        # {results}[result_array]
    
    for Corr, Op in zip(Correlators, Op_names):
        result["C_"+Op] = Corr
    return result

def calc_corr_tilde(Correlators, args):                               # Correlators [Ops][N_T] 
    Op_names = args[0]                                                       # Initialization of args   # same length as Correlators
    if len(Correlators) != len(Op_names):
        print("length of Corr is not same as Ops, calc_corr_tilde")
        return 1
    result = {}                                                 # {results}[result_array]
    
    for Corr, Op in zip(Correlators, Op_names):
        result["C_tilde_"+Op] = []
        for i in range(len(Corr)-2):                                        # 2 oder 1, I have to decide
            result["C_tilde_"+Op].append(Corr[i]-Corr[i+2])
    return result

def calc_eff_mass_log(Correlators, args):                               # Correlators [Ops][N_T] 
    Op_names = args[0]                                                       # Initialization of args   # same length as Correlators
    if len(Correlators) != len(Op_names):
        print("length of Corr is not same as Ops, calc_eff_mass_log")
        return 1
    result = {}                                                 # {results}[result_array]
    
    for Corr, Op in zip(Correlators, Op_names):
        result["m_eff_log_"+Op] = []
        for i in range(len(Corr)-1):
            result["m_eff_log_"+Op].append(np.log(Corr[i])/np.log(Corr[i+1]))
    return result

# def calc_eff_mass_impl(Correlators, args):                               # Correlators [Ops][N_T] 
#     Op_names = args[0]                                                       # Initialization of args   # same length as Correlators
#     if len(Correlators) != len(Op_names):
#         print("length of Corr is not same as Ops, calc_eff_mass_impl")
#         return 1
#     result = {}  

#     def zero_eff_mass(eff_mass, ratio, index):
#         return np.cosh(eff_mass*(T_2-index))/np.cosh(eff_mass*(T_2-(index+1))) - ratio                                               # {results}[result_array]
    
#     for Corr, Op in zip(Correlators, Op_names):
#         result["m_eff_impl_"+Op] = []
#         for i in range(len(Corr)-1):                                        # 2 oder 1, I have to decide
#             ratio = Corr[i]/Corr[i+1]
#             T_2 = (len(Corr)-2)//2
#             result["m_eff_impl_"+Op].append(bisect(f=zero_eff_mass, a=1e-30, b=100, args = (ratio,i)))
#     return result

def calc_eff_mass_impl_deri(Correlators, args):                               # Correlators [Ops][N_T] 
    Op_names = args[0]                                                       # Initialization of args   # same length as Correlators
    if len(Correlators) != len(Op_names):
        print("length of Corr is not same as Ops, calc_eff_mass_impl_deri")
        return 1
    result = {}  

    def zero_eff_mass(eff_mass, ratio, index):
        return np.sinh(eff_mass*(T_2-index))/np.sinh(eff_mass*(T_2-(index+1))) - ratio                                               # {results}[result_array]
    for Corr, Op in zip(Correlators, Op_names):
        result["m_eff_impl_deri_"+Op] = []
        for i in range(len(Corr)-3):                                        # 2 oder 1, I have to decide
            ratio = (Corr[i]-Corr[i+2])/(Corr[i+1]-Corr[i+3])
            T_2 = (len(Corr)-2)//2
            if (T_2-i) == 0 or (T_2-(i+1)) == 0:                    # Only happens for sinh. For cosh both values are well-defined
                # result["m_eff_impl_deri_"+Op].append(float("inf"))
                result["m_eff_impl_deri_"+Op].append(0)
            else:
                res = bisect(f=zero_eff_mass, a=1e-30, b=1000, args = (ratio,i))
                if np.isnan(res):
                    result["m_eff_impl_deri_"+Op].append(0)
                else:
                    result["m_eff_impl_deri_"+Op].append(bisect(f=zero_eff_mass, a=1e-30, b=1000, args = (ratio,i)))
    return result
    
def calc_convexity(Correlators, args):                          # Correlators [Ops][N_T] 
    Op_names = args[0]                                                       # Initialization of args   # same length as Correlators
    if len(Correlators) != len(Op_names):
        print(len(Correlators),len(Op_names))
        print("length of Corr is not same as Ops, calc_corr")
        return 1
    result = {}                                        # {results}[result_array]
    
    for Corr, Op in zip(Correlators, Op_names):
        result["convexity_"+Op] = []
        for i in range(len(Corr)-2):
            result["convexity_"+Op].append((Corr[i]-2*Corr[i+1]+Corr[i+2])/4)
    return result


def basic_analysis(Correlators, args):                                                      # Give C_pi, C_rho, C_pipi
    result = {}
    names = args[0]
    if len(Correlators)%3 != 0:
        print("Len of Corr in basic_analysis is not a mutiple of 3. ")
        return 1
    if len(names) != len(Correlators)//3:
        print("Len of Correlators in basic_analysis is not len of names*3. ")
        return 2
    Op_names = ["pi", "rho", "pipi"]
    for key, value in calc_corr(Correlators, [Op_names,]).items():
        result[key]=value
    for key, value in calc_corr_tilde(Correlators, [Op_names,]).items():
        result[key]=value
    for key, value in calc_eff_mass_log(Correlators, [Op_names,]).items():
        result[key]=value
    for key, value in calc_eff_mass_impl_deri(Correlators, [Op_names,]).items():
        result[key]=value
    for key, value in calc_convexity(Correlators, [Op_names,]).items():
        result[key]=value
    # for key, value in calc_eff_mass_impl(Correlators, [Op_names,]).items():
    #     result[key]=value
    for key in result:
        print(key, result[key])
    return result