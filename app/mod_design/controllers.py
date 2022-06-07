import pandas as pd
import numpy as np
import os
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.backends.backend_pdf
#import matplotlib.pyplot as plt
#import seaborn as sns

from app.mod_design.models import runAutomation

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(APP_ROOT, '../resources/')

def validateFile(filename):

	return '.' in filename and filename.rsplit('.', 1)[1].lower() in set(['xlsx', 'xls', 'csv'])

def readFile(path_to_fluo, path_to_od):

	fluos = pd.read_csv(path_to_fluo, index_col='index')
	ods = pd.read_csv(path_to_od, index_col='index')

	return fluos, ods

def automateDesign(t, libraries, od_params):

	result = runAutomation(t, libraries, od_params)

	return "Hello, world!"

	'''
	params = [par['Value'] for par in libraries]

	
	for ar in result[:5]:

		circuit = Circuit(fluo_selected, od_selected, gate_params_selected, od_params_selected, None, ar[0])
		model = gate_wrapper2(f_t, circuit)
		gfp = pd.Series(model[33], index=f_t/60)



		plt.plot(gfp, label=ar[0])
		mar_solutions.append((ar[0], gfp))
	plt.legend(bbox_to_anchor=(1.05, 1), ncol=2)
	sns.despine()
	plt.show()

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
	'''
