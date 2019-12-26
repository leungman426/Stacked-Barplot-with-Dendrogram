# input file :
# count of heterozygous and homozygous SNPs (count_snp)
# SNP profiling for all samples(new_sample_df)
# SNPs genotype table for case (case_df)

# output: show the heterzoygous/homozggous SNPs level for all the samples
# and the hierarchical clustering of the samples' SNP genotype in the same plot

count_snp = pd.read_csv('/count_snp', sep=" ", stringsAsFactors = FALSE)

# hierarchical clustering matrix
tran_new_sample_df = new_sample_df.transpose()
y = pdist(tran_new_sample_df, metric = "cityblock", p = 2)
z = linkage(y, method = 'complete')
index_phylogenetics = dendrogram(z)['ivl']

# sort the samples in the same way as the phylogenetic plot
def sortdf(df,**kwargs):
    new_df = pd.DataFrame()
    if 'colname' in kwargs.keys():
        k1 = kwargs['colname']
        for item in k1:
            selection = df[df['samples'] == item]
            new_df = new_df.append(selection)

    if 'index' in kwargs.keys():
        k2 = kwargs['index']
        for item in k2:
            selection = df[df.index == int(item)]
            new_df = new_df.append(selection)
    new_df.index = range(0,len(new_df))
    return new_df

plot_df = sortdf(count_snpe, index=index_phylogenetics)

# stacked barplot for heterozygous and homozygous SNPs count
ypos = np.arange(len(plot_df['samples']))
plt.figure(figsize=(70, 30), dpi = 90)
plt.subplots_adjust(hspace=0.6, right=1, left=0.05)
plt.subplot(211)
a1 = plt.bar(ypos, plot_df['heterozygous'], color="b", bottom=plot_df['homozygous'])
a2 = plt.bar(ypos, plot_df['homozygous'], color="g")
plt.xticks([], [])
plt.ylim(0,75)
plt.xlim(0,len(ypos))
plt.legend([a1,a2], ['heterozygous', 'homozygous'], fontsize=50)
plt.ylabel('Number of the Risk SNPs', fontsize=50)

# label the sample name in red for the dendrogram
def labelcolor(labelnames,col):
    ax = plt.gca()
    for item in ax.get_xticklabels():
        for item1 in labelnames:
            if item.get_text() == item1:
                item.set_color(col)

# hierarchical clustering
plt.subplot(212)
ax = dendrogram(z, orientation="bottom", labels=plot_df['samples'].tolist(), leaf_rotation=90, leaf_font_size=36)
dendrogram(z, orientation="bottom", distance_sort=True, labels=plot_df['samples'].tolist(), leaf_rotation=90, leaf_font_size=36)
labelnames = list(case_df)
labelcolor(labelnames, 'r')
