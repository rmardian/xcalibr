import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

import pymc3 as pm
from pymc3.ode import DifferentialEquation
import arviz as az
import theano
THEANO_FLAGS = "optimizer=fast_compile"


class SGNModel:
    
    @staticmethod
    def hill_equation(x, K, n):
        return x**n / (K**n + x**n)
    
    @staticmethod
    def gfp_only_model(y, t, p):
        
        #dependent variables
        Auto, OD = y[0], y[1]
        #a = p[0]
        #alpha, beta = extra
        a, alpha, beta = p[0], p[1], p[2]
        
        gamma = SGNModel.growth_rate(t, OD, alpha, beta)
        #differential equations
        dOD = gamma * OD
        dAuto = a - gamma * Auto
        return [dAuto, dOD]
    
    @staticmethod
    def gate_model_no_auto(y, t, p):
        
        bn, bc, bg, syn_ECFn, syn_ECFc, syn_ECF, deg, syn_GFP, deg_GFP, K, n = p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10]
        ECFn, ECFc, ECF, GFP = y[0], y[1], y[2], y[3]
        ind1, ind2 = 1, 1
        
        #gamma = growth_rate(t, OD, alpha, beta)
        #differential equations
        #dOD = gamma * OD
        dECFn = bn + syn_ECFn * ind1 - deg * ECFn
        dECFc = bc + syn_ECFc * ind2 - deg * ECFc
        dECF = syn_ECF * ECFn * ECFc - deg * ECF
        
        dGFP = bg + syn_GFP * SGNModel.hill_equation(ECF, K, n) - deg_GFP * GFP

        return [dECFn, dECFc, dECF, dGFP]

fluos = pd.read_csv('marionette_fluo.csv', index_col='time')
ods = pd.read_csv('marionette_od.csv', index_col='time')
gates = list(set([i[:-3] for i in fluos.columns.tolist()]))

gate = 'e11x32STPhoRadA'

fluo_sel = fluos.loc[:, fluos.columns.str.startswith(gate)]
od_sel = ods.loc[:, ods.columns.str.startswith(gate)]

fluo = fluo_sel.iloc[:,3]
od = od_sel.iloc[:,3]


pars = {
    'bn': 9.980960e+00,
    'bc': 4.873768e+00,
    'bg': 5.160116e+00,
    'syn_ECFn': 9.989766e+01,
    'syn_ECFc': 5.087405e+01,
    'syn_ECF': 9.844930e-09,
    'syn_GFP': 9.351119e-03,
    'deg': 3.627036e+03,
    'deg_GFP': 1.347098e-01,
    'K': 4.908509e+01,
    'n': 2.490464e+00
}
errs = {
    'bn': 1.282754e-01,
    'bc': 1.242285e+01,
    'bg': 1.022647e+00,
    'syn_ECFn': 4.283182e+00,
    'syn_ECFc': 2.649231e+01,
    'syn_ECF': 1.152230e-08,
    'syn_GFP': 9.366010e-03,
    'deg': 9.459877e-03,
    'deg_GFP': 2.046354e-01,
    'K': 3.432406e+01,
    'n': 2.512038e+00
}

with pm.Model() as od_model:
    
    bn = pm.Normal('bn', mu=pars['bn'], sigma=errs['bn'])
    bc = pm.Normal('bc', mu=pars['bc'], sigma=errs['bc'])
    bg = pm.Normal('bg', mu=pars['bg'], sigma=errs['bg'])
    syn_ECFn = pm.Normal('syn_ECFn', mu=pars['syn_ECFn'], sigma=errs['syn_ECFn'])
    syn_ECFc = pm.Normal('syn_ECFc', mu=pars['syn_ECFc'], sigma=errs['syn_ECFc'])
    syn_ECF = pm.Normal('syn_ECF', mu=pars['syn_ECF'], sigma=errs['syn_ECF'])
    syn_GFP = pm.Normal('syn_GFP', mu=pars['syn_GFP'], sigma=errs['syn_GFP'])
    deg = pm.Normal('deg', mu=pars['deg'], sigma=errs['deg'])
    deg_GFP = pm.Normal('deg_GFP', mu=pars['deg_GFP'], sigma=errs['deg_GFP'])
    K = pm.Normal('K', mu=pars['K'], sigma=errs['K'])
    n = pm.Normal('n', mu=pars['n'], sigma=errs['n'])
    
    y_hat = DifferentialEquation(
        func=SGNModel.gate_model_no_auto, times=fluo.index, n_states=4, n_theta=11
    )(y0=[0, 0, 0, 0], theta=[bn, bc, bg, syn_ECFn, syn_ECFc, syn_ECF, deg, syn_GFP, deg_GFP, K, n])
    
    fluo_est = pm.Normal('fluo', mu=y_hat.T[0], sd=0.25, observed=fluo)
    
    #trace = pm.sample(2000, tune=1000, cores=1)
    
    step = pm.Metropolis()
    #trace = pm.sample(2000, tune=1000, cores=1, step=step)
    trace = pm.sample(1000, tune=1000, cores=1, step=step)