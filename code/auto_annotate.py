import scanpy as sc
import numpy as np
import pandas as pd
from tqdm import tqdm
sc.settings.verbosity = 0
import warnings
warnings.filterwarnings('ignore')

def auto_annot(adata, cluster, db = 'Panglao',
               species = 'Mm', cutoff = 0.25):
    # Read anndata object
    data = adata.copy() 

    if db == 'Panglao':
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

    if db == 'CellMarker':

        if (species == 'Mm'):
            # Load CellMarker Database
            db = pd.read_csv('https://github.com/Xenon8778/Auto_cell_annot/raw/main/data/Cell_marker_Human.csv')
            db_sub = db[['speciesType','geneSymbol','cellName','tissueType']]

            # Selecting species and subsetting database
            db_sub = db_sub[(db_sub['species'] == species) | (db_sub['species'] == 'Mm Hs') | (db_sub['species'] == 'Hs Mm')]
            db_sub = db_sub[db_sub['ubiquitousness index'] < 0.1]

            db_prep = []
            for i in np.unique(db_sub['cell_name']):
                x1 = [x for x in db_sub[db_sub['cell_name'] == i]['Symbol']]
                x2 = [i,x1]
                db_prep.append(x2)
            
        else:
            db = pd.read_csv('https://github.com/Xenon8778/Auto_cell_annot/raw/main/data/Cell_marker_Human.csv')

            # Selecting species and subsetting database
            db_sub = db_sub[(db_sub['species'] == species) | (db_sub['species'] == 'Mm Hs') | (db_sub['species'] == 'Hs Mm')]
            db_sub = db_sub[db_sub['ubiquitousness index'] < 0.1]

            db_prep = []
            for i in np.unique(db_sub['cell_name']):
                x1 = [x for x in db_sub[db_sub['cell_name'] == i]['Symbol']]
                x2 = [i,x1]
                db_prep.append(x2)

    

    # Annotate each cluster 
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

    # df_out

    # Selected the cell type with max score
    maxValueIndex = df_out.idxmax(axis=1)
    res = pd.concat([maxValueIndex,df_out.max(axis=1)], axis = 1)
    res.columns = ['Cell Type','Score']
    res['Cell Type Conf'] = [res['Cell Type'][i] if res['Score'][i]> cutoff else 'Unknown' for i in res.index]
    return(res)