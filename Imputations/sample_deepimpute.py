import scanpy 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import logging
import scprep
from deepimpute.multinet import MultiNet


logging.basicConfig(filename='log_impute.txt',level=logging.INFO)


#Loading dataset
logging.info('Loading dataset')
adata = scanpy.read_h5ad("/home/je334775/Impute/Spleen_gene_count.h5ad")


#converting to dataframe
df = adata.to_df()


#implementing deepimute
logging.info('Implementing deepimpute')
multinet = MultiNet() 
# Using custom parameters
NN_params = {
        'learning_rate': 1e-4,
        'batch_size': 64,
        'max_epochs': 200,
        'ncores': 5,
        'sub_outputdim': 512,
            'architecture': [
            {"type": "dense", "activation": "relu", "neurons": 200},
            {"type": "dropout", "activation": "dropout", "rate": 0.3}]
    }

multinet = MultiNet(**NN_params)


#Using all of the data for training
logging.info('Training the model')
multinet.fit(df)


logging.info('Performing Imputation on the dataset')
imputedData = multinet.predict(df)


df_zeros = df[df==0.0].count().sum()
print('Number of zeros in the original dataset:',df_zeros)


impute_zeros = imputedData[imputedData==0.0].count().sum()
print('Number of zeros in the Imputed dataset:',impute_zeros)
print('Percent imputed: ',((df_zeros-impute_zeros)/df_zeros)*100)

logging.info('Visualizing')
limits = [0,100]
df = df[:500]
imputedData = imputedData[:500]

fig,ax = plt.subplots()

jitter = np.random.normal(0,1,df.size) # Add some jittering to better see the point density
ax.scatter(df.values.flatten()+jitter,imputedData.values.flatten(),s=2)
ax.plot(limits,limits,'r-.',linewidth=2)
ax.set_xlim(limits)
ax.set_ylim(limits)

plt.show()
plt.savefig('/home/je334775/Impute/fig1.png', bbox_inches='tight', dpi=300)


multinet.test_metrics


fig, (ax1, ax2) = plt.subplots(1,2, figsize=(16,6))

columns_name_orig = df.columns
columns_name_imp = imputedData.columns

scprep.plot.scatter(x=df[columns_name_orig[0]], y=df[columns_name_orig[1]], c=df[columns_name_orig[2]],
                    ax=ax1, xlabel='first gene', ylabel='second gene', legend_title='third gene', title='Before DeepImpute')
scprep.plot.scatter(x=imputedData[columns_name_orig[0]], y=imputedData[columns_name_orig[1]], c=imputedData[columns_name_orig[2]],
                    ax=ax2, xlabel='first gene', ylabel='second gene', legend_title='third gene', title='After DeepImpute')
plt.show()
plt.tight_layout()
plt.savefig('/home/je334775/Impute/fig2.png', bbox_inches='tight', dpi=300)







