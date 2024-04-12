# scAutoCT
Automatic annotation of single-cells post-clustering in Python using the PanglaoDB and CellMarker2.0 database.

## Intructions
Pull the git folder into your working directory 
```bash
git clone https://github.com/Xenon8778/scAutoCT
```

Loading and running the function
```python
from scAutoCT.code.auto_annotate import auto_annot
res  = auto_annot(data, cluster='leiden', species='Mm')
new_cluster_names = list(res['Cell Type Conf'])
```
Replacing Leiden (or any other clustering algorithm labels) with cell types -
```python

Cell_dict = dict(zip(list(np.arange(0, len(np.unique(data.obs['leiden'])),1)), new_cluster_names))
data.obs['Auto_labels'] = data.obs['leiden']
data.obs['Auto_labels'] = data.obs['Auto_labels'].map(Cell_dict)
res
```
## Databases
- PangaloDB. https://doi.org/10.1093/database/baz046
- CellMarker2.0. https://doi.org/10.1093/nar/gkac947
