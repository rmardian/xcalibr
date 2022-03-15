import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.integrate import odeint
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

par = {
    'bn': 1 * 10**-1,
    'bc': 1 * 10**-1,
    'bg': 1 * 10**-2,
    'syn_ECFn': 4 * 10**0,
    'syn_ECFc': 4 * 10**0,
    'syn_ECF': 50 * 10**-10,
    'deg': 7 * 10**-3,
    'syn_GFP': 1 * 10**4,
    'deg_GFP': 1 * 10**-2,
    'a': 1 * 10**-2,
    'K': 1 * 10**-1,
    'n': 2 * 10**0
}

############### STATIC MODELS ###############

def growth_rate(t, OD, alpha, beta):
    return (alpha * (1 - (OD/beta)))

def od_wrapper(t, r, c, c0):
    
    def od_model(OD, t, r, c):
        dOD = growth_rate(t, OD[0], r, c) * OD[0]
        return dOD
    
    od_rates = (r, c)
    od_sol = odeint(od_model, c0, t, od_rates)
    return od_sol[:,0]

def od_inference(od):
    od_bounds = [(0, 0, 0), (1, 2, 0.1)]
    od_params = []
    od_t = od.index
    for idx in range(od.shape[1]):
        od_data = od.iloc[:,idx]
        opt, _ = curve_fit(od_wrapper, od_t, od_data, bounds=od_bounds)
        od_params.append(opt)
    return od_params

def hill_equation(x, K, n):
    return x**n / (K**n + x**n)

############### ODE MODELS ###############

#model with auto-fluorescence term
def gate_wrapper_complete(t, bn, bc, bg, syn_ECFn, syn_ECFc, syn_ECF, deg, syn_GFP, deg_GFP, a, K, n, ind1, ind2, extra, y0):
    
    #fixed parameters
    alpha, beta = extra

    def gate_model(y, t):
        
        #dependent variables
        ECFn, ECFc, ECF, GFP, Auto, OD = y
        
        gamma = growth_rate(t, OD, alpha, beta)
        #differential equations
        dOD = gamma * OD
        dECFn = bn + syn_ECFn * ind1 - (deg + gamma) * ECFn
        dECFc = bc + syn_ECFc * ind2 - (deg + gamma) * ECFc
        dECF = syn_ECF * ECFn * ECFc - (deg + gamma) * ECF
        
        dGFP = bg + syn_GFP * hill_equation(ECF, K, n) - (deg_GFP + gamma) * GFP
        dAuto = a - gamma * Auto

        return [dECFn, dECFc, dECF, dGFP, dAuto, dOD]
    
    solution = odeint(gate_model, y0, t)
    return solution.transpose()

############### SAMPLING INITIAL GUESSES ###############

#objective function
def computeSSE(init_params, f_t, f_data, bounds, od_params, ind1, ind2):
    
    extra = (od_params[0], od_params[1])
    y0 = pd.Series(np.append(np.zeros(5), od_params[2]))
    
    def model_fit(t, a, b, c, d, e, f, g, h, i, j, k, l):
        fit = gate_wrapper_complete(t, a, b, c, d, e, f, g, h, i, j, k, l, ind1, ind2, extra, y0)
        return pd.Series(fit[3]) + pd.Series(fit[4])
    
    f_params, _ = curve_fit(model_fit, f_t, f_data, p0=init_params, bounds=bounds)
    solution = model_fit(f_t, *f_params)
    error = [(val-sal)**2 for val, sal in zip(solution, f_data)]
    
    return sum(error)

#generate random numbers from a uniform distribution for the initial guesses
def randomSearch(iterations, num_params, f_t, f_data, bounds, od_params, ind1, ind2):
    
    initialGuesses = []
    for k in tqdm(range(iterations)):
        guess = [np.random.uniform(low=low, high=high) for low, high in zip(bounds[0], bounds[1])]
        error = computeSSE(guess, f_t, f_data, bounds, od_params, ind1, ind2)
        initialGuesses.append((error, guess))
    
    sortedGuesses = sorted(initialGuesses)
    bestGuess = sortedGuesses[0][1]
    
    return bestGuess

############### INFERENCE MODULES ###############

def fit_single_state(fluo, od_params, ind1, ind2, optimize=False):
    
    parameters = list(par.keys())
    f_t = fluo.index
    f_data = fluo.copy()
    extra = (od_params[0], od_params[1])
    y0 = pd.Series(np.append(np.zeros(5), od_params[2]))
    
    if optimize:
        #with sampling initial guesses
        lower_bounds = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        upper_bounds = [1e1, 1e1, 1e1, 1e2, 1e2, 1e-3, 1e-0, 1e5, 1e0, 1e1, 1e2, 5]
        gate_bounds = [lower_bounds, upper_bounds]
        init_guesses = randomSearch(10, len(parameters), f_t, f_data, gate_bounds, od_params, ind1, ind2)
    else:
        #without sampling initial guesses
        init_guesses = [par[i] for i in parameters]
        lower_bounds = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        upper_bounds = [1e1, 1e1, 1e1, 1e2, 1e2, 1e-4, 1e-1, 1e5, 1e0, 1e1, 1e2, 4]
        gate_bounds = [lower_bounds, upper_bounds]
    
    def model_fit(t, a, b, c, d, e, f, g, h, i, j, k, l):
        fit = gate_wrapper_complete(t, a, b, c, d, e, f, g, h, i, j, k, l, ind1, ind2, extra, y0)
        return pd.Series(fit[3]) + pd.Series(fit[4])
    
    f_params, f_cov = curve_fit(model_fit, f_t, f_data, p0=init_guesses, bounds=gate_bounds)
    f_sim = model_fit(f_t, *f_params)
    f_df = pd.DataFrame({'Parameters': parameters, 'Value': f_params, 'Err': np.sqrt(np.diag(f_cov))})
    return f_df, f_sim, f_t, f_data

def fit_all_states(fluo, od_params):
    
    num_states = 4
    num_vars = 5
    
    title = fluo.columns.tolist()[0] #take the first entry only because they are all the same
    f_t = np.concatenate([fluo.index] * num_states)
    f_data = pd.concat([fluo.iloc[:,i] for i in range(num_states)]) 
    y0 = pd.concat([pd.Series(np.append(np.zeros(num_vars), od_params[i][2])) for i in range(num_states)])
    
    #gate_bounds = [lower_bounds, upper_bounds]
    #init_guesses = randomSearch_global(10, len(parameters), f_t, f_data, gate_bounds, od_params)
    
    parameters = list(par.keys())
    init_guesses = list(par.values())
    lower_bounds = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    upper_bounds = [1e1, 1e1, 1e1, 1e2, 1e2, 1e-4, 1e-1, 1e5, 1e0, 1e1, 1e2, 4]
    gate_bounds = [lower_bounds, upper_bounds]
    
    def model_fit_inner(t, a, b, c, d, e, f, g, h, i, j, k, l, ind1, ind2, extra, y0):
        fit = gate_wrapper_complete(t, a, b, c, d, e, f, g, h, i, j, k, l, ind1, ind2, extra, y0)
        return pd.Series(fit[3]) + pd.Series(fit[4])
    
    def model_fit(t, a, b, c, d, e, f, g, h, i, j, k, l):
        result = [model_fit_inner(fluo.index, a, b, c, d, e, f, g, h, i, j, k, l, int(n/2), n%2, \
                                   (od_params[n][0], od_params[n][1]), np.append(np.zeros(num_vars), od_params[n][2])) \
                                   for n in range(num_states)]
        return pd.concat(result)
    
    f_params, f_cov = curve_fit(model_fit, f_t, f_data, p0=init_guesses, bounds=gate_bounds)#, method='dogbox')
    #f_params, f_cov = curve_fit(model_fit, f_t, f_data, p0=init_guesses)
    f_sim = model_fit(f_t, *f_params)
    f_df = pd.DataFrame({'Parameters': parameters, 'Value': f_params, 'Err': np.sqrt(np.diag(f_cov))})
    return f_df, f_sim, f_t, f_data

def fit_all_gates(fluos, ods, gates, all_states=False):

	f_dfs, f_sims, f_ts, f_datas = [], [], [], []

	for sel, gate in tqdm(enumerate(gates)):
		
		gate = gates[sel]
		fluo_sel = fluos.loc[:, fluos.columns.str.startswith(gate)]
		od_sel = ods.loc[:, ods.columns.str.startswith(gate)]
		
		od_params_sel = od_inference(od_sel)
		
		if not all_states:
			n = 3
			f_df, f_sim, f_t, f_data = fit_single_state(fluo_sel.iloc[:,n], od_params_sel[n], int(n/2), n%2) #only 11 state
		if all_states:
			f_df, f_sim, f_t, f_data = fit_single_state(fluo_sel.iloc[:,n], od_params_sel[n])
		
		f_dfs.append(f_df)
		f_sims.append(f_sim)
		f_ts.append(f_t)
		f_datas.append(f_data)

	return f_dfs, f_sims, f_ts, f_datas