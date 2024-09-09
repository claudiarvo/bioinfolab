# scRNASeq Analysis (bioinfolab)

## Overview
This project is designed for the analysis of single-cell RNA sequencing (scRNA-seq) data, with support for various imputation methods and clustering approaches. It provides functionalities such as preprocessing, dimensionality reduction, imputation, clustering, and evaluation using metrics like Adjusted Rand Index (ARI) and silhouette scores.

## Key Features
- Multiple Imputation Methods: Choose from methods like MAGIC, DeepImpute, KNN, and SCVI for data imputation.
- Dimensionality Reduction: Supports PCA and t-SNE for visualizing high-dimensional data.
- Clustering: Implements Leiden clustering with adjustable resolution parameters.
- Evaluation: Computes ARI and silhouette scores to assess the quality of the clustering.

## Pipeline Folder

#### workflow.py
This is the main script to run the scRNA-seq pipeline. It takes in arguments for the data path, imputation method, dataset name, and clustering resolution. 

Key Arguments:

* --data_path: Path to the input .h5ad file containing scRNA-seq data.
* --imputation: Imputation method to use. Options: MAGIC, DeepImpute, KNN, SCVI.
* --dataset_name: Name of the dataset (e.g., Thymus, Spleen).
* --resolution: Clustering resolution (default: 0.1).

```
python3 workflow.py --data_path /path/to/dataset.h5ad --imputation DeepImpute --dataset_name Thymus --resolution 0.01

```

#### scRNAAnalysis.py

This script contains the scRNAAnalysis class, which performs the various steps in the scRNA-seq analysis workflow, such as loading data, preprocessing, imputation, dimensionality reduction, clustering, and evaluation.

Key Methods:

* load_data(): Loads the dataset from the given file path.
* preprocess_data(): Filters and normalizes the data.
* imputation(): Applies the selected imputation method.
* run_dimension_reduction(): Reduces the dataset dimensions using PCA and t-SNE.
* clustering(): Clusters the data using the Leiden algorithm and saves clustering plots.
* evaluation(): Evaluates clustering quality using silhouette scores and ARI.

## Imputations Folder

This folder contains scripts that implement various imputation methods for scRNA-seq data. These methods handle missing or zero-inflated values in the gene expression matrix, and the scripts save both the imputed data and visualizations comparing the data before and after imputation.


#### autoimpute_data.py
This script applies the AutoImpute to perform imputation on scRNA-seq data. The imputation strategy (e.g., mean, median, mode) can be passed as an argument.

Arguments:

* --strategy: Imputation strategy (mean, median, most_frequent, etc.)

Workflow:

* Reads the data from a .h5ad file and transposes it.
* Applies the specified imputation strategy using SingleImputer.
* Saves the imputed data to a CSV file.

```
python3 autoimpute_data.py --strategy mean

```

#### magic_imputation.py

This script uses the MAGIC algorithm to impute missing values in scRNA-seq data. It reads the input data from an .h5ad file and then applies MAGIC to a subset of the data.

Workflow:

* Reads the data from an .h5ad file.
* Subsets the data for faster processing.
* Applies the MAGIC algorithm to the dataset.
* Saves the imputed data and generates a plot comparing the data before and after imputation.

```
python3 magic_imputation.py

```

#### sample_deepimpute.py
This script applies the DeepImpute method to perform deep learning-based imputation on scRNA-seq data. DeepImpute handles the zero-inflated nature of scRNA-seq datasets and predicts missing values using a neural network model.

Workflow:

* Loads the dataset from an .h5ad file.
* Configures and trains a neural network model (MultiNet) with custom parameters.
* Imputes missing values from the dataset.
* Compares the number of zeros before and after imputation.
* Visualizes the imputation results using scatter plots and saves the results.

```
python3 sample_deepimpute.py

```

## Dependencies

To run the pipeline, you need the following Python libraries:

* scanpy
* magic-impute
* numpy
* pandas
* matplotlib
* scikit-learn
* deepimpute
* scvi-tools

Install the required packages using the following command:

```
pip install -r requirements.txt

```

## How to Use

- Prepare your data: Ensure your scRNA-seq data is in the .h5ad format.
- Run the pipeline: Execute the workflow.py script with appropriate arguments.
- Check outputs: Logs and clustering plots will be saved to the pipeline_output/logs and pipeline_output/leiden_* directories, respectively.

## Output 

- Detailed logs of the pipeline execution will be saved in the pipeline_output/logs/ directory.
- Clustering results before and after imputation are saved as PNG images in the pipeline_output/leiden_before_imputation/ and pipeline_output/leiden_after_imputation/ directories.