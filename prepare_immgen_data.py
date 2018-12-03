## downloaded from immgen.org, go to processed data request

# GEO gives details on normalization, but no guarantees processed data on immgen web site
# were processed the same way:

#The datasets were pre-filtered to keep only those probesets for which a gene symbol 
#could be found in the Affymetrix annotation. CEL files were normalized using Affymetrix 
#Power Tools on the predefined probeset ID list mentioned above, and using the standard RMA workflow
#(background adjustment, quantile normalization, median polish probeset summarization). Output on 2^ scale.

# the processed data download page on immgen does NOT give details on normalization
# but I don't think it's log!

def quantileNormalize(df_input,ref=None):
    """
    df_input - genes (or probes) x samples
    
    """
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    if type(ref)!=type(None):
        rank = sorted_df[ref].tolist()
    else:
        rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df

path1 = 'immgen_microarray_downloaded_May2018/RequestedImmGenData2018-05-14_21-35-07.csv'
path2 = 'immgen_microarray_downloaded_May2018/RequestedImmGenData2018-05-14_21-36-50.csv'

phase1 = pd.read_csv(path1)
phase2 = pd.read_csv(path2)

phase1.columns = [i.strip(' ').strip('`') for i in phase1.columns]
phase2.columns = [i.strip(' ').strip('`') for i in phase2.columns]
                  
# confirm that first 3 columns are the same (product over all booleans)
print np.product(phase1.iloc[:,:3].values == phase2.iloc[:,:3].values)
    
# concatenate
immgen = pd.concat([phase1,phase2.iloc[:,3:]],axis=1)
print phase1.shape, phase2.shape, immgen.shape

# all probe ids are unique:
print len(immgen['ProbeSetID'].unique()) == immgen.shape[0]

# therefore can be safely moved to index:
immgen.index = immgen['ProbeSetID']
immgen.drop('ProbeSetID',inplace=True,axis=1)

# Are description and genesymbol 1-to-1? No, but there are fewer description, just drop them
print len(immgen['GeneSymbol'].unique()),len(immgen['Description'].unique())
immgen.drop('Description',inplace=True,axis=1)


#slow below

#collapse probes by max mean as described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3166942/pdf/1471-2105-12-322.pdf
print 'collapsing...'
# mask for selecting duplicates
isdup = immgen.duplicated('GeneSymbol',keep=False).values

# unique separately
immgen[~isdup].shape

dupgen = immgen[isdup].copy()
groups = dupgen.groupby('GeneSymbol').groups

collapsed = {}
counter = 0
start = time.time()
for key,value in groups.items():
    counter+=1
    # get the probe id with max mean expression
    maxmean = dupgen.drop('GeneSymbol',axis=1).loc[value].mean(axis=1).idxmax()
    collapsed[key] = dupgen.drop('GeneSymbol',axis=1).loc[maxmean]
    if counter/500.==counter/500:
        print counter,'/',len(groups)
        print (time.time()-start)/60.,'min.'
print (time.time()-start)/60.,'min.'
collapsed = pd.DataFrame(collapsed)

# select probes that were unique
uq = immgen[~isdup].copy()

uq.index = uq['GeneSymbol']
uq.drop('GeneSymbol',axis=1,inplace=True)

# concatenate and overwrite "immgen"
cat = pd.concat([uq,collapsed.T])

# order index alphabetically
cat = cat.loc[sorted(cat.index)]

# quantile normalize. I decide to normalize after excluding probes I am not even considering
print 'quantile normalizing...'
cat = quantileNormalize(cat)

# turn into tpms + 10
cat = (cat*1e6/cat.sum()) + 10

# white spaces in gene names!
cat.index = [i.strip(' ') for i in cat.index]

cat.to_csv('backups/immgen_max_mean_collapsed_quant_norm_tpm_10.tsv',sep='\t')