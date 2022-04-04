import pandas as pd
import numpy as np
import gseapy
import pickle


def get_bp_dataset():
	bp_datasets = ["E-GEOD-46449", "E-GEOD-5388", "E-GEOD-5392", "E-GEOD-53987", "GSE12649"]

	for dataset in bp_datasets:
		data = pd.read_table(f"../data/bipolar/{dataset}_data.tsv", index_col=0).transpose().reset_index()

		temp_cols = ['index']
		temp_cols.extend([col.replace('"', '').split('.')[0] for col in data.columns[1:]])
		data.columns = temp_cols

		des = pd.read_table(f"../data/bipolar/{dataset}_design.tsv", index_col=0)

		des_map = {
			des.at[row, 'FileName'].replace('"', '').split('.')[0]: des.at[row, 'Target']
			for row in des.index
		}

		cls = [des_map[pat] for pat in data.columns[1:]]

		yield (dataset, data, cls)


def get_sc_dataset():
	sc_datasets = ['E-GEOD-12649', 'E-GEOD-21138', 'E-GEOD-21935', 'E-GEOD-53987', 'GSE93987']

	for dataset in sc_datasets:
		data = pd.read_table(f"../data/Scz/{dataset}_data.tsv", index_col=0).transpose().reset_index()

		temp_cols = ['index']
		temp_cols.extend([col.replace('"', '').replace('.CEL', '') for col in data.columns[1:]])
		data.columns = temp_cols

		des = pd.read_table(f"../data/Scz/{dataset}_design.tsv", index_col=0)

		des_map = {
			des.at[row, 'FileName'].replace('"', '').replace('.CEL', ''): des.at[row, 'Target'].replace('schizophrenic', 'schizophrenia')
			for row in des.index
		}

		cls = [des_map[pat] for pat in data.columns[1:]]

		yield (dataset, data, cls)
        

def get_dm_dataset():
	sc_datasets = ["GSE13760", "GSE15653", "GSE16415", "GSE20966", "GSE23343"]

	for dataset in sc_datasets:
		data = pd.read_table(f"../data/T2DM/{dataset}_data.tsv", index_col=0).transpose().reset_index()

		temp_cols = ['index']
		temp_cols.extend([col.replace('"', '').replace('.CEL', '') for col in data.columns[1:]])
		data.columns = temp_cols

		des = pd.read_table(f"../data/T2DM/{dataset}_design.tsv", index_col=0)

		des_map = {
			des.at[row, 'FileName'].replace('"', '').split('.')[0]: des.at[row, 'Target']
			for row in des.index
		}

		cls = [des_map[pat] for pat in data.columns[1:]]

		yield (dataset, data, cls)


def main():
	for tup in get_sc_dataset():
		dataset, data, cls = tup

		print(f"Processing Schizophrenia Dataset: {dataset}")

		gsea = gseapy.gsea(
			data=data,
			gene_sets="kegg.gmt",
			cls=cls,
			max_size=500,
			min_size=15,
			permutation_num=100,
			outdir=None,
			no_plot=True,
			processes=1,
			format='png'
		)

		with open(f'gsea_res/sc_{dataset}.pkl', 'wb') as f:
			pickle.dump(gsea, f)


	for tup in get_bp_dataset():
		dataset, data, cls = tup

		print(f"Processing Bipolar Dataset: {dataset}")

		gsea = gseapy.gsea(
			data=data,
			gene_sets="kegg.gmt",
			cls=cls,
			max_size=500,
			min_size=15,
			permutation_num=100,
			outdir=None,
			no_plot=True,
			processes=1,
			format='png'
		)

		with open(f'gsea_res/bp_{dataset}.pkl', 'wb') as f:
			pickle.dump(gsea, f)
            
	for tup in get_dm_dataset():
		dataset, data, cls = tup

		print(f"Processing T2DM Dataset: {dataset}")

		gsea = gseapy.gsea(
			data=data,
			gene_sets="kegg.gmt",
			cls=cls,
			max_size=500,
			min_size=15,
			permutation_num=100,
			outdir=None,
			no_plot=True,
			processes=1,
			format='png'
		)

		with open(f'gsea_res/dm_{dataset}.pkl', 'wb') as f:
			pickle.dump(gsea, f)

if __name__ == '__main__':
	main()
