from autoimpute.imputations import SingleImputer
import pandas as pd
import scanpy as sc
import logging
import argparse

# Read arguments
parser = argparse.ArgumentParser(description="Strategy for SimpleImputer")
parser.add_argument('--strategy', type=str, help='Enter the strategy')
args = parser.parse_args()

# Read and transpose data 
logging.info('Loading dataset...')
adata_gex = sc.read_h5ad("Spleen_gene_count.h5ad")
logging.info('Transposing dataset...')
data = adata_gex.to_df().T

# Example using default instance of MiceImputer
logging.info('Implementing SingleImputer...')
if args.strategy != None:
    si = SingleImputer(strategy=args.strategy) # pass through data onceg
else:
    si = SingleImputer(strategy='mean') # Mean strategy by default

# fit transform returns a generator by default, calculating each imputation method lazily
logging.info('Training the model and transforming data...')
data_imputed = si.fit_transform(data[:10])

logging.info('Saving data imputed: ')
data_imputed.to_csv('data_imputed.csv')
#print(data_imputed)
