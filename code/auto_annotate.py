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

    # Selecting Database

    ## PanglaoDB database PMID: 30951143
    if db == 'Panglao':
        print('Using PanglaoDB')
        # Load PanglaoDB Database
        db = pd.read_table('https://raw.githubusercontent.com/Xenon8778/Auto_cell_annot/main/data/PanglaoDB_markers_27_Mar_2020.tsv', sep = '\t')
        db_sub = db[['species','official gene symbol','cell type','ubiquitousness index']]

        # Selecting species and subsetting database
        db_sub = db_sub[(db_sub['species'] == species) | (db_sub['species'] == 'Mm Hs') | (db_sub['species'] == 'Hs Mm')]
        db_sub = db_sub[db_sub['ubiquitousness index'] < 0.1]

        db_prep = []
        for i in np.unique(db_sub['cell type']):
            x1 = [x for x in db_sub[db_sub['cell type'] == i]['official gene symbol']]
            if species == 'Mm':
                x1 = [x.capitalize() for x in x1]
            x2 = [i,x1]
            db_prep.append(x2)
        print('Database Prepped')

    ## CellMarker2.0 database PMID: 36300619
    elif db == 'CellMarker':
        print('Using CellMarker 2.0')
        
        # Mouse 
        if (species == 'Mm'):
            # Load CellMarker Database
            db = pd.read_csv('https://github.com/Xenon8778/Auto_cell_annot/raw/main/data/Cell_marker_Mouse.csv')
            db_sub = db[['species','Symbol','cell_name','tissue_class']]
            print(db_sub['tissue_class'].unique())

            if tissuelist == 'All': 
                db_sub = db_sub
            else:
                db_sub = db_sub[db_sub['tissue_class'].isin(tissuelist)]

            db_prep = []
            for i in np.unique(db_sub['cell_name']):
                x1 = [x for x in db_sub[db_sub['cell_name'] == i]['Symbol']]
                x1 = [j for j in x1 if str(j) != 'nan']
                if len(x1) > 1:
                    x2 = [i,x1]
                    db_prep.append(x2)

        # Human or others    
        else:
            db = pd.read_csv('https://github.com/Xenon8778/Auto_cell_annot/raw/main/data/Cell_marker_Human.csv')
            db_sub = db[['species','Symbol','cell_name','tissue_class']]
            print(db_sub['tissue_class'].unique())

            if tissuelist == 'All': 
                db_sub = db_sub
            else:
                db_sub = db_sub[db_sub['tissue_class'].isin(tissuelist)]

            db_prep = []
            for i in np.unique(db_sub['cell_name']):
                x1 = [x for x in db_sub[db_sub['cell_name'] == i]['Symbol']]
                x1 = [j for j in x1 if str(j) != 'nan']
                if len(x1) > 1:
                    x2 = [i,x1]
                    db_prep.append(x2)
        print('Database Prepped')

    # Annotating each cluster 
    for i in tqdm(range(len(db_prep))):
        try:
            sc.tl.score_genes(data, db_prep[i][1], score_name= db_prep[i][0])
        except:
            pass
        try:
            df = data.obs[[cluster,db_prep[i][0]]]
            if i == 0:
                df_out = df.groupby(cluster).mean()
            else:
                df_out = pd.concat([df_out,df.groupby(cluster).mean()], axis= 1)
        except:
            pass

    # Selected the cell type with max score
    maxValueIndex = df_out.idxmax(axis=1)
    res = pd.concat([maxValueIndex,df_out.max(axis=1)], axis = 1)
    res.columns = ['Cell Type','Score']
    res['Cell Type Conf'] = [res['Cell Type'][i] if res['Score'][i]> cutoff else 'Unknown' for i in res.index]
    return(res)