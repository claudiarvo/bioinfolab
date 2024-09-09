import numpy as np
import pandas as pd
import scanpy as sc
import scprep
import magic
import matplotlib.pyplot as plt
import logging
logging.basicConfig(filename='log_trial.txt',level=logging.INFO)

#read data
anndata = sc.read_h5ad('/home/hg621061/biolab/project/data/Spleen_gene_count.h5ad')
logging.info(f'Anndata: {anndata}')
#transform to df
df = anndata.to_df()
df = df[:500]
#df.to_csv('/home/hg621061/biolab/project/data/Spleen_gene_count.csv')
logging.info('Original data saved as csv')
#init magic object and apply to the data
magic_op = magic.MAGIC()
df_magic = magic_op.fit_transform(df)
#df_magic.to_csv('/home/hg621061/biolab/project/data/imputed_magic_Spleen_gene_count.csv')
logging.info('Imputed data saved as csv')

logging.info('Visualizing')
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(16,6))

columns_name_orig = df.columns
columns_name_imp = df_magic.columns

scprep.plot.scatter(x=df[columns_name_orig[0]], y=df[columns_name_orig[1]], c=df[columns_name_orig[2]],
                    ax=ax1, xlabel='first gene', ylabel='second gene', legend_title='thirds gene', title='Before MAGIC')
scprep.plot.scatter(x=df_magic[columns_name_orig[0]], y=df_magic[columns_name_orig[1]], c=df_magic[columns_name_orig[2]],
                    ax=ax2, xlabel='first gene', ylabel='second gene', legend_title='thirds gene', title='After MAGIC')
plt.tight_layout()
plt.savefig('/home/hg621061/biolab/project/data/magic_imputation.png', bbox_inches='tight', dpi=300)
