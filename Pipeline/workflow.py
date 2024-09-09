# COMMAND: python3 workflow.py --data_path /work/lect0094/RNA/h5ad/Thymus.h5ad --imputation DeepImpute --dataset_name Thymus --resolution 0.01

import logging
import argparse
from scRNAAnalysis import scRNAAnalysis

# Parse arguments
parser = argparse.ArgumentParser(description="scRNA analysis pipeline")
parser.add_argument('--data_path', type=str, help='Enter the data path', required=True)
parser.add_argument('--imputation', type=str, choices=['MAGIC', 'DeepImpute', 'KNN', 'SCVI'], default='MAGIC', help='Enter the method you want to use for imputation.')
parser.add_argument('--dataset_name', type=str, help='Enter the name of the dataset, e.g. Spleen', required=True)
parser.add_argument('--resolution', type=float, help='Choose the resolution ', default=0.1, required=True)
args = parser.parse_args()

# Prepare log file
logging.basicConfig(filename = f'/home/lect0094/group3/pipeline_output/logs/workflow_{args.imputation}_{args.dataset_name}.log', filemode='w',
                    level = logging.INFO,
                    format = '%(asctime)s:%(levelname)s:%(name)s:%(message)s')

# Create an instance of the class and run analysis
analysis = scRNAAnalysis(args.data_path, args.imputation, args.dataset_name, args.resolution)
analysis.run_analysis()
