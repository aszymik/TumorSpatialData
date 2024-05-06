import pandas as pd 
import matplotlib.pyplot as plt 
from scipy.cluster.hierarchy import dendrogram, linkage


# Df preprocessing
columns = ['patient', 'component_size', 'other', 'CD15-Tumor', 'CD15+Tumor', 'Tcell', 'Bcell', 'BnTcell', 'Macrophage', 'Neutrophil', 'DC']
df = pd.read_csv('if_data/cell_type_percentages_in_TLSs.tsv', sep='\t', index_col=0)
df = df.reindex(columns=columns)
df = df.drop('component_size', axis=1)

# Hclust
linked = linkage(df.drop('patient', axis=1), 'ward')

# Plot
fig, axes = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={'width_ratios': [1, 2]})

# Dendrogram
dendro = dendrogram(linked,
                    orientation='left',
                    distance_sort='descending',
                    no_labels=True,
                    show_leaf_counts=False,
                    ax=axes[0])

sample_order = dendro['ivl']
df.index = df.index.astype(str)
df = df.reindex(sample_order)

# Bar plot
df.plot( 
	x = 'patient', 
	kind = 'barh', 
	stacked = True,
    ax=axes[1]
)

# Adjust axes
for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.yaxis.tick_right()
    ax.tick_params(bottom=False, right=False)

axes[0].set_ylabel('')
axes[0].tick_params(labelbottom=False)

axes[1].set_facecolor('#f1f1f1')
axes[1].set_xlabel('Cell percentage')
axes[1].set_xlim(0, 1)
axes[1].set_yticklabels(ax.get_yticklabels(), fontsize=8)
axes[1].set_ylabel('Patient')
axes[1].yaxis.set_label_position("right")

handles, labels = axes[1].get_legend_handles_labels()
axes[1].legend().remove()
fig.legend(handles, labels, loc='upper left', fontsize=7, bbox_to_anchor=(0.015, 0.92))

fig.suptitle('Cell distribution in TLS candidates')
plt.subplots_adjust(wspace=0.01)
plt.show()
