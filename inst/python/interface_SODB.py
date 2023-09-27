import pysodb
import anndata
import squidpy

def get_SODB_dataset(dataset_name: str,
                     experiment_name: str = "default" ) -> anndata.AnnData:

    '''
    INPUT: 
        dataset_name - name of dataset, must exist within the SODB
        experiment_name - name of experiment, must exist within the given dataset. See details.
    OUTPUT:
        adata - anndata object populated with the desired data
    
    DETAILS:
        
        Default behavior is to use the first associated experiment name with a provided dataset
        if no experiment_name is provided. 
        
        If a provided dataset or experiment is not found, an empty anndata object will be returned.
        This will cause downstream function calls to list out available datasets and
        experiments.
    
    
    '''
    
    sodb = pysodb.SODB()
    
    if dataset_name not in sodb.list_dataset():
        adata = anndata.AnnData()
        adata.uns["dataset_not_found"] = False
        return (adata)
    
    if experiment_name == "default":
        experiment_name = sodb.list_experiment_by_dataset(dataset_name)[0]
    elif experiment_name not in sodb.list_experiment_by_dataset(dataset_name):
        adata = anndata.AnnData()
        adata.uns["experiment_not_found"] = False
        return (adata)

    
    adata = sodb.load_experiment(dataset_name, experiment_name)

    return adata

def check_SODB_adata(dataset_name: str,
                     adata: anndata.AnnData,
                     experiment_name: str = "default"):
    '''
    INPUT: 
        dataset_name - name of dataset provided to upstream function `get_SODB_dataset()`
        experiment_name - name of experiment provided to upstream function `get_SODB_dataset()`
        adata - anndata object outoput from upstream function `get_SODB_dataset()`
    OUTPUT:
        NONE: checks for errors in anndata object
    DETAILS:
        
        If a provided dataset or experiment is not found, this function will throw an error
        which will list the available datasets or the available experiments for a particular dataset.
        
    
    
    '''
    uns_keys = list(adata.uns.keys())
    
    if ("dataset_not_found" not in uns_keys and "experiment_not_found" not in uns_keys):
        pass
    else:
        sodb = pysodb.SODB()
        if("dataset_not_found" in uns_keys):
            
            avail_datasets = sodb.list_dataset()
            
            err_msg = f"Provided dataset name '{dataset_name}' not found. Using Giotto, run `listSODBDatasetNames()` to view all available names."
            
            raise KeyError(err_msg)
        elif("experiment_not_found" in uns_keys):
            avail_experiments = sodb.list_experiment_by_dataset(dataset_name)
            err_msg = f"Provided experiment name '{experiment_name}' not found within dataset '{dataset_name}'. The available experiements for this dataset are as follows:{avail_experiments} "
            raise KeyError(err_msg)

def list_SODB_datasets(category:str) -> list:
    avail_names = ["Spatial Transcriptomics", 
                   "Spatial Proteomics",
                   "Spatial Metabolomics",
                   "Spatial Genomics",
                   "Spatial MultiOmics",
                   "All"]
        
    if category not in avail_names:
        err_msg = f"Provided category '{category}' not found within the available names, which are as follows:{avail_names}."
        raise KeyError(err_msg)
    
    sodb = pysodb.SODB()
    
    if category == "All": 
        return sodb.list_dataset()
    else:
        return sodb.list_dataset_by_category(category)
    
def list_SODB_dataset_experiments(dataset_name:str) -> list:
    
    sodb = pysodb.SODB()

    avail_datasets = sodb.list_dataset()
    
    if dataset_name not in avail_datasets:
        err_msg = f"Provided dataset name '{dataset_name}' not found. Using Giotto, run `listSODBDatasetNames()` to view all available names."
        raise KeyError(err_msg)
    
    return sodb.list_experiment_by_dataset(dataset_name)
