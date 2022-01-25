
# import warnings
# warnings.simplefilter(action='ignore', category=OutdatedPackageWarning)

import numpy as np
import pandas as pd
import pingouin as pg
from statsmodels.stats import multitest

from pathos import multiprocessing



def _get_mediation_analysis_results(data): 

	mediation_results = {}
	r2_results = {}

	for i,row in data.iterrows(): 
		df = row.unstack(0)  # dataframe with sample rows and omic columns
		has_atac = not df["atac"].isna().all()
		has_erna = not df["erna"].isna().all()
		has_rna = not df["rna"].isna().all()

		try: 
			if has_atac and has_erna: 
				mediators = ["atac", "erna"]
			elif has_atac: 
				mediators = ["atac"]
			elif has_erna: 
				mediators = ["erna"]
			else: 
				pass

			mediation_results[i] = pg.mediation_analysis(data=df, x="snps", y="rna", m=mediators)
		except: 
			pass

		_r2_results = {}
		try:
			_r2_results["atac_erna_med"] = pg.partial_corr(data=df, x="snps", y="rna", covar=["atac", "erna"]).iloc[0]
		except:
			pass
		try:
			_r2_results["atac_med"] = pg.partial_corr(data=df, x="snps", y="rna", covar=["atac"]).iloc[0]
		except:
			pass
		try:
			_r2_results["erna_med"] = pg.partial_corr(data=df, x="snps", y="rna", covar=["erna"]).iloc[0]
		except:
			pass
		try:
			_r2_results["direct"] = pg.partial_corr(data=df, x="snps", y="rna").iloc[0]
		except:
			pass
		if len(_r2_results) > 0: 
			r2_results[i] = pd.concat(_r2_results)


		# try:
		# 	_r2_results = {}
		# 	if has_atac and has_erna: 
		# 		_r2_results["atac_erna_med"] = pg.partial_corr(data=df, x="snps", y="rna", covar=["atac", "erna"]).iloc[0]
		# 	if has_atac: 
		# 		_r2_results["atac_med"] = pg.partial_corr(data=df, x="snps", y="rna", covar=["atac"]).iloc[0]
		# 	if has_erna: 
		# 		_r2_results["erna_med"] = pg.partial_corr(data=df, x="snps", y="rna", covar=["erna"]).iloc[0]

		# 	_r2_results["direct"] = pg.partial_corr(data=df, x="snps", y="rna").iloc[0]

		# 	r2_results[i] = pd.concat(_r2_results)
		# except: 
		# 	pass

	if len(mediation_results) > 0:
		mediation_results = pd.concat(mediation_results)
		mediation_results = mediation_results.droplevel(1).pivot(columns="path").swaplevel(0,1,axis=1).sort_index(axis=1)
	else: 
		mediation_results = pd.DataFrame()

	if len(r2_results) > 0:
		r2_results = pd.concat(r2_results).unstack([1,2])
	else: 
		r2_results = pd.DataFrame()

	return mediation_results, r2_results


def get_mediation_analysis_results(data): 

	split_data = np.array_split(data, 16)
	with multiprocessing.ProcessingPool(16) as pool: 
		results = pool.map(_get_mediation_analysis_results, split_data)

	mediation_results, r2_results = zip(*results)
	mediation_results = pd.concat(mediation_results)
	r2_results = pd.concat(r2_results)

	# return mediation_results, r2_results

	# Fix na values when only one mediator is tested
	if "Indirect" in mediation_results.columns.levels[0]:
		for col in mediation_results.columns.levels[1]: 
			mediation_results.loc[~mediation_results["Y ~ atac", col].isna(), ("Indirect atac", col)] = mediation_results.loc[~mediation_results["Y ~ atac", col].isna(), ("Indirect atac", col)].fillna(mediation_results["Indirect", col])
			mediation_results.loc[~mediation_results["Y ~ erna", col].isna(), ("Indirect erna", col)] = mediation_results.loc[~mediation_results["Y ~ erna", col].isna(), ("Indirect erna", col)].fillna(mediation_results["Indirect", col])
		mediation_results.drop(columns="Indirect", level=0, axis=1, inplace=True)
	mediation_results.columns = mediation_results.columns.remove_unused_levels()

	# return mediation_results, r2_results

	return mediation_results, r2_results


