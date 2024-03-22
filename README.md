# Auto_cell_annot
Automatic annotation of single-cells after clustering in python using Scanpy and PanglaoDB

## Steps to run
```python
from auto_annotate import auto_annot
res  = auto_annot(data, cluster='leiden', species='Mm')
new_cluster_names = list(res['Cell Type'])
```
Replacing leiden (or any other clustering algorithm labels) with cell types -
```python

Cell_dict = dict(zip(list(np.arange(0, 34,1)), new_cluster_names))
data.obs['Auto_labels'] = data.obs['leiden']
data.obs['Auto_labels'] = data.obs['Auto_labels'].astype(int).map(Cell_dict)
res
```
