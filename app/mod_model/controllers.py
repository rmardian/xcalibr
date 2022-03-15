import pandas as pd
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import seaborn as sns

from app.mod_model.models import fit_all_gates

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(APP_ROOT, '../resources/')

def validateFile(filename):

	return '.' in filename and filename.rsplit('.', 1)[1].lower() in set(['xlsx', 'xls', 'csv'])

def readFile(path_to_fluo, path_to_od):

	fluos = pd.read_csv(path_to_fluo, index_col='index')
	ods = pd.read_csv(path_to_od, index_col='index')

	return fluos, ods

def generateModel(fluos, ods, gates):

	f_dfs, f_sims, f_ts, f_datas = fit_all_gates(fluos, ods, gates, all_states=False)
	return f_dfs, f_sims, f_ts, f_datas

def plot_fitting(f_ts, f_datas, f_sims, gates):

	f, axs = plt.subplots(2, 6, sharex=True, sharey=True, figsize=(14, 4))
	axr = axs.ravel()
	for i, ax in enumerate(axr):
		if i < len(gates):
			ax.scatter(f_ts[i]/60, f_datas[i], c='slategrey', s=7)
			ax.plot(f_ts[i]/60, f_sims[i], c='deeppink')
			ax.set_title(gates[i])
			ax.set_xlabel('Time (h)')
		else:
			ax.set_visible(False)
	sns.despine()
	plt.suptitle("Curve Fitting - 11 State")
	plt.tight_layout()

	return f

def plot_parameters(f_df, gates, show_err=False):
    
    parameters = f_df[0]['Parameters'].values
    n_paras = len(parameters)
    values = np.stack([df['Value'].values for df in f_df])
    errors = np.stack([df['Err'].values for df in f_df])
    
    f, axs = plt.subplots(int(n_paras/3), 3, sharex=True, figsize=(16, (n_paras+1)/2))
    for i, ax in enumerate(axs.ravel()):
        if i < n_paras:
            if show_err:
                ax.errorbar(gates, values[:,i], errors[:,i], fmt='o', label='experiment')
            else:
                ax.scatter(gates, values[:,i], label='experiment')
            ax.set_ylabel(parameters[i])
        else:
            ax.set_visible(False)
        #ax.set_xticks(rotation=90) 
        ax.set_xticklabels(gates, rotation=90)
    sns.despine()
    plt.suptitle("Parameters Analysis")
    plt.tight_layout()

def generatePlot(f_ts, f_datas, f_dfs, f_sims, gates):

	pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join(RESOURCES, "plot.pdf"))
	f = plot_fitting(f_ts, f_datas, f_sims, gates)
	pdf.savefig(f)

	f = plot_parameters(f_dfs, gates)
	pdf.savefig(f)

	f = plot_parameters(f_dfs, gates, True)
	pdf.savefig(f)

	pdf.close()

def generateReport(f_dfs, gates):

	df = pd.DataFrame()
	for i, f_df in enumerate(f_dfs):
		f_df['Gate'] = gates[i]
		df = df.append(f_df)

	df[['Gate', 'Parameters', 'Value', 'Err']].to_csv(os.path.join(RESOURCES, "report.csv"), index=False)
