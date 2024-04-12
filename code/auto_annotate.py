import scanpy as sc
import numpy as np
import pandas as pd
from tqdm import tqdm
sc.settings.verbosity = 0
import warnings
warnings.filterwarnings('ignore')

def auto_annot(adata, cluster, db = 'Panglao', tissuelist = ['All'],
               species = 'Hs', cutoff = 0.25):
    
    # data -> input annadata object with cluster labels
    # cluster -> anndata 'obs' column containing the cluster labels
    # db -> Database to be used
    # species -> Specify species to use
    # tissuelist -> Input tissue types from CellMarker 2.0 Dataset in a list
    # cutoff -> Score threshold for confidence of annotation


    # Read anndata object
    data = adata.copy() 

    # Scale data
    print('Scaling...')
    sc.pp.scale(data)

    # Selecting Database

    ## PanglaoDB database PMID: 30951143
    if db == 'Panglao':
        if species == 'Mm':
            print('Using PanglaoDB')
            file = pd.read_csv('../data/PanglaoDB_2020_Mouse.csv')
            
            #Subsetting by Tissue
            if 'All' in tissuelist:
                file = file
            else:
                file = file[file['Tissue'].isin(tissuelist)].reset_index()

            db_prep = [] 
            for i in range(file.shape[0]):
                a = file['CT'][i]
                b = file['Gene'][i].split(',') 
                db_prep.append([a,b])
            db_prep
            print('Database Prepped')
        else: 
            print('Using PanglaoDB')
            file = pd.read_csv('../data/PanglaoDB_2020_Human.csv')
            
            #Subsetting by Tissue
            if 'All' in tissuelist:
                file = file
            else:
                file = file[file['Tissue'].isin(tissuelist)].reset_index()
            
            db_prep = [] 
            for i in range(file.shape[0]):
                a = file['CT'][i]
                b = file['Gene'][i].split(',') 
                db_prep.append([a,b])
            db_prep
            print('Database Prepped')

    ## CellMarker2.0 database PMID: 36300619
    elif db == 'CellMarker':
        print('Using CellMarker 2.0')
        
        # Mouse 
        if (species == 'Mm'):
            file = pd.read_csv('../data/CellMarker2_Mouse_Aug.csv')
            
            #Subsetting by Tissue
            if 'All' in tissuelist:
                file = file
            else:
                file = file[file['Tissue'].isin(tissuelist)].reset_index()

            db_prep = [] 
            for i in range(file.shape[0]):
                a = file['CT'][i]
                b = file['Gene'][i].split(',') 
                db_prep.append([a,b])
            db_prep

        # Human or others    
        else:
            file = pd.read_csv('../data/CellMarker2_Human_Aug.csv')
            
            #Subsetting by Tissue
            if 'All' in tissuelist:
                file = file
            else:
                file = file[file['Tissue'].isin(tissuelist)].reset_index()
            
            db_prep = [] 
            for i in range(file.shape[0]):
                a = file['CT'][i]
                b = file['Gene'][i].split(',') 
                db_prep.append([a,b])
            db_prep
            
        print('Database Prepped')

    # Annotating each cluster 
    df_out = pd.DataFrame(index= adata.obs[cluster].unique())
    for i in tqdm(range(len(db_prep))):
        try:
            sc.tl.score_genes(data, db_prep[i][1], score_name= db_prep[i][0], ctrl_size = 100,
                                n_bins = 24, random_state = 1)
        except:
            pass
        try:
            df = data.obs[[cluster,db_prep[i][0]]]
            df_out[db_prep[i][0]] = df.groupby(cluster).mean()
        except:
            pass

    # Selected the cell type with max score
    maxValueIndex = df_out.idxmax(axis=1)
    res = pd.concat([maxValueIndex,df_out.max(axis=1)], axis = 1)
    res.columns = ['Cell Type','Score']
    res['Cell Type Conf'] = [res['Cell Type'][i] if res['Score'][i]> cutoff else 'Unknown' for i in res.index]
    return(res)