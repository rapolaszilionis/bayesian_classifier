def bayesian_classifier(op,cp):
    '''
    op - observed gene expression profile, genes x samples
    cp - class profiles, genes x samples, same genes as op
    returns log10(P(E|type)), the max value is the closes cell type
    '''
    
    #we assume that each cell type has a well define centroid, let's represent this expression vector
    #as the fractions of all mRNAs for each genes (i.e. normalized the expression such that the expression of
    #all genes sums to 1)
    
    cp = cp/cp.sum()
    
    #we assume that the exact expression pattern we observe (E) is multinomially distributed around the centroid.
    #Bayes' formula: P(type|E) = P(E|type)*P(type)/P(E)
    #our classifier is naive, so each E is equally likely (this is how I interpret "naive", although it
    #may have more to do with the assumption that genes are uncorrelated)
    
    ptes = pd.DataFrame({cell:(np.log10(cp.T.values)*op[cell].values).sum(axis=1) for cell in op.columns})
    ptes.index = cp.columns
    return ptes
    
#by Adrian Veres for saving and loading 
def save_df(obj, filename):
    np.savez_compressed(filename, data=obj.values, index=obj.index.values, columns=obj.columns.values)
    
def load_df(filename,encoding=u'ASCII'):
	u"""you may want to specify encoding='latin1'
	when loading python 2 pickle with python 3.
	https://stackoverflow.com/questions/28218466/unpickling-a-python-2-object-with-python-3
	"""
	with np.load(filename,encoding=encoding) as f:
		obj = pd.DataFrame(**f)
	return obj
	
##########################################################################################
 
# cat - pd.DataFrame with immgen data   
cat = pd.read_csv('backups/immgen_max_mean_collapsed_quant_norm_tpm_10.tsv',sep='\t',index_col=0)
cat.head()

#get a reduced dense dataframe for this purpose
# gene_list - your gene_list
common_genes = [i for i in cat.index if i in gene_list]
print len(common_genes)

mask_genes = np.array([i in common_genes for i in gene_list])
print mask_genes.sum()

start = time.time()
bays = []
i = 0
for j in range(5000,Eraw.shape[0]+5000,5000):
    
    # Eraw - sparse cells x gene matrixx
    j = min(j,Eraw.shape[0])
    tmp_dense = pd.DataFrame(Eraw.T[mask_genes][:,i:j].todense())
    tmp_dense.index = np.array(gene_list)[mask_genes]
    
    bay = bayesian_classifier(tmp_dense,cat.loc[common_genes])
    bays.append(bay)
    i = j
    
    print(time.time()-start)/60.,'min.'
    print 'cells from %d to %d done'%(i-5000,min(j,Eraw.shape[0]))
    
bay = pd.concat(bays,axis=1)
save_df(bay,'backups/bayesian_immgen_results')

#19 min to run on 15k cells