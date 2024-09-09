import scanpy as sc
import numpy as np
import magic
import pickle
import pandas as pd
import logging
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score, adjusted_rand_score
from deepimpute.multinet import MultiNet
from knn_smoothing import knn_smoothing
import scvi

class scRNAAnalysis:

    def __init__(self, data_path, imputation_method, dataset_name, resolution):
        self.data_path = data_path
        self.imputation_method = imputation_method
        self.dataset_name = dataset_name
        self.adata = None
        self.adata_df = None
        self.resolution = resolution

    def load_data(self):
        logging.info(f'Loading {self.dataset_name} dataset...')
        self.adata = sc.read_h5ad(self.data_path)
        logging.info(f"Raw adata: {self.adata}")

    def preprocess_data(self):
        
        # Filtering, normalization, log transformation, etc.
        logging.info('Starting preprocess...')

        sc.pp.calculate_qc_metrics(self.adata, percent_top=None, log1p=False, inplace=True)
        
        # Filtering cells
        sc.pp.filter_cells(self.adata, min_genes=self.adata.obs['n_genes_by_counts'].min()+1) 
        # Filtering genes
        sc.pp.filter_genes(self.adata, min_cells=int(self.adata.var['n_cells_by_counts'].min()) + 1)

        # Filtering outliers
        q1 = np.percentile(self.adata.obs.total_counts, 25)
        q3 = np.percentile(self.adata.obs.total_counts, 75)
        cutoff_n_genes = q3 + 1.5*(q3-q1)
        logging.info(f'Cut off value n genes: {cutoff_n_genes}')
        self.adata = self.adata[self.adata.obs.total_counts < cutoff_n_genes, :]
        self.adata_df = self.adata.to_df() #DeepImpute takes raw data as input
        # Normalization
        sc.pp.normalize_total(self.adata)
        # Log
        sc.pp.log1p(self.adata)
        # highly variables
        sc.pp.highly_variable_genes(self.adata)

        logging.info('Preprocess completed')

    def imputation(self):
        
        if self.imputation_method=='SCVI':
            logging.info('Imputation using SCVI...')
            scvi.model.SCVI.setup_anndata(self.adata)
            model = scvi.model.SCVI(self.adata)
            model.to_device('cuda:0') # for using GPU in training
            model.train()
            data_imputed = model.get_normalized_expression()
            adata_imputed = sc.AnnData(data_imputed, data_imputed.index.to_frame(), data_imputed.columns.to_frame())
            logging.info('Imputation completed')
            return adata_imputed

        if self.imputation_method=='MAGIC':
            logging.info('Imputation using MAGIC...')
            imputer = magic.MAGIC()
            adata_imputed = imputer.fit_transform(self.adata)
            logging.info('Imputation completed')
            return adata_imputed
        
        if self.imputation_method=='DeepImpute': 
            logging.info('Imputation using Deepimpute...')
            multinet = MultiNet()
            logging.info('Training the model and transforming data...')
            subset = 1
            VMR = self.adata.var['dispersions'].min()
            NN_params = {
                    
                    'learning_rate': 1e-4,
                    'batch_size': 128,
                    'max_epochs': 100,                    
                    'sub_outputdim': 1500,
            }
            multinet = MultiNet(**NN_params)
            multinet.fit(self.adata_df,cell_subset = subset, minVMR = VMR)
            data_imputed = multinet.predict(self.adata_df)
            logging.info('Imputation completed...')
            logging.info('Transforming imputed df to Anndata Object...')
            adata_imputed = sc.AnnData(data_imputed, data_imputed.index.to_frame(), data_imputed.columns.to_frame())
            logging.info('Imputation completed')
            logging.info('Normalizing and logarithmizing imputed data for dim reduction and clustering')
            # Normalization
            sc.pp.normalize_total(adata_imputed)
            # Log
            sc.pp.log1p(adata_imputed)
            # highly variables
            sc.pp.highly_variable_genes(adata_imputed)
            return adata_imputed

    def run_dimension_reduction(self, adata):

        # Perform dimensionality reduction
        logging.info('Starting imensionality Reduction...')
        sc.tl.pca(adata, n_comps=100)
        sc.tl.tsne(adata, n_pcs = 50)
        logging.info('Dimensionality Reduction completed')

    def clustering(self, adata, save_path, title, case, resolution):

        sc.pp.neighbors(adata, use_rep = 'X_tsne')
        logging.info(f'NEIGHBORS done')

        logging.info(f'Clustering with LEIDEN {resolution}...')
        logging.info('')
        sc.tl.leiden(adata, resolution = resolution, key_added = f'leiden_{case}') #default resolution in 1.0

        sc.pl.tsne(adata, color=[f'leiden_{case}'], title=title, legend_loc = 'on data')
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        logging.info('clustering figure saved')
        logging.info('')

    def compute_ARI(self, adata, case):

        # Compute ARI
        logging.info('Computing ARI...')

        if case=='baseline':
            ari = adjusted_rand_score(self.adata.obs.Main_cluster_name[self.adata.obs.Main_cluster_name != 'NA'], self.adata.obs.leiden_baseline[self.adata.obs.Main_cluster_name != 'NA'])
            logging.info(f'Leiden labels: {self.adata.obs.leiden_baseline.values}')
        else:
            ari = adjusted_rand_score(self.adata.obs.Main_cluster_name[self.adata.obs.Main_cluster_name != 'NA'], adata.obs.leiden_model[self.adata.obs.Main_cluster_name != 'NA'])
            logging.info(f'Leiden labels: {adata.obs.leiden_model.values}')

        logging.info(f'Model ARI score: {round(ari,4)}')

    def evaluation(self, adata, case):
        logging.info('Starting evaluation...')
        
        if case=='baseline':
            silhouette = silhouette_score(adata.obsm['X_tsne'], adata.obs.leiden_baseline.values)
            logging.info(f'Baseline imputation: Silhouette leiden: {round(silhouette,4)}')
            
        else:
            silhouette = silhouette_score(adata.obsm['X_tsne'], adata.obs.leiden_model.values)
            logging.info(f'Model imputation: Silhouette leiden: {round(silhouette,4)}')

    def baseline_analysis(self):
        self.preprocess_data()
        self.run_dimension_reduction(self.adata)
        save_path = f"/home/lect0094/group3/pipeline_output/leiden_before_imputation/leiden_baseline_{self.dataset_name}.png"
        self.clustering(self.adata, save_path, f'{self.dataset_name} clustering Baseline', case='baseline', resolution=self.resolution)
        self.evaluation(self.adata, case='baseline')
        self.compute_ARI(self.adata, case="baseline")
    

    def model_analysis(self):
        adata_imputed = self.imputation()
        self.run_dimension_reduction(adata_imputed)
        save_path = f"/home/lect0094/group3/pipeline_output/leiden_after_imputation/leiden_model_imputation_{self.imputation_method}_{self.dataset_name}.png"
        self.clustering(adata_imputed, save_path, f'{self.dataset_name} clustering  Model',  case='model', resolution=self.resolution)
        self.evaluation(adata_imputed, case='model')
        self.compute_ARI(adata_imputed, case='model')

    def run_analysis(self):

        self.load_data()
        self.baseline_analysis()
        self.model_analysis()
