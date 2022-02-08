
# coding: utf-8

# In[1]:


import Bio.SeqIO, sys, pandas, re, argparse, os


# In[2]:


if sys.argv[1] == "-f":
    mode = "debug"
else:
    mode = "batch"
    
if mode=='debug':
    args = {
        'file_cds':'/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/cds/Sus_scrofa.Sscrofa11.1.cds.all.fa',
        'file_id':'/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/id_mapping/sscrofa_gene_ensembl.tsv',
        'file_out':'/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/cds/Sus_scrofa.Sscrofa11.1.cds.all.out.fa',
        'translate_table':1,
        'remove_stop':1,
    }
    os.chdir('/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/cds')
elif mode=='batch':
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_cds', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--file_id', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--file_out', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--translate_table', metavar='INT', default=0, type=int, help='NCBI translate table number. 0 for no translation.')
    parser.add_argument('--remove_stop', metavar='INT', default=0, type=int, help='')
    args = parser.parse_args()  
    g = dict()
    for attr in [a for a in dir(args) if not a.startswith('_')]:
        g[attr] = getattr(args, attr)
    args = g


# In[3]:


df_id = pandas.read_csv(args['file_id'], sep='\t')
df_id = df_id.loc[:,['Gene stable ID', 'Transcript stable ID']].drop_duplicates()
df_id.columns = ['gene_id','transcript_id']

seqids = list()
seqlens = list()
with open(args['file_cds'], "r") as handle:
    for record in Bio.SeqIO.parse(handle, "fasta"):
        seqids.append(record.id)
        seqlens.append(len(record.seq))
df = pandas.DataFrame([seqids, seqlens]).T
print('# transcripts in input fasta:', df.shape[0])

df.columns = ['transcript_id','length']
df['transcript_id'] = df['transcript_id'].str.replace('\\..*', '')
df = pandas.merge(df, df_id, on='transcript_id', how='left')
print('# genes with multiple transcripts:', df.loc[(df['gene_id'].duplicated()),'gene_id'].drop_duplicates().shape[0])
max_length = pandas.DataFrame(df.groupby('gene_id')['length'].max())
max_length['gene_id'] = max_length.index
df = pandas.merge(max_length, df, how='left')
is_surplus_longest = (df['gene_id'].duplicated())
print('# genes with multiple longest transcripts:', df.loc[is_surplus_longest,'gene_id'].drop_duplicates().shape[0])
df = df.loc[(~is_surplus_longest),:]
print('# longest transcripts:', df.shape[0])


# In[4]:


transcript_ids = df['transcript_id'].tolist()
seqs = list()
not_multiple_of_3 = list()
with open(args['file_cds'], "r") as handle:
    for record in Bio.SeqIO.parse(handle, "fasta"):
        transcript_id = re.sub('\\..*', '', record.id)
        if (transcript_id in transcript_ids):
            gene_id = df.loc[df['transcript_id']==transcript_id,'gene_id'].values[0]
            record.id = gene_id
            record.description = ''
            record.name = ''
            if args['translate_table']!=0:
                is_multiple_of_3 = (len(record.seq)%3==0)
                seq_protein = record.seq.translate(table=args['translate_table'])
                seq_protein._data = re.sub('\\*$', '', seq_protein._data)
                if not is_multiple_of_3:
                    not_multiple_of_3.append(record.id)
                    #print('sequence not multiple of 3:', record.id, seq_protein)
                if args['remove_stop']:
                    seq_protein._data = seq_protein._data.replace('*', '')
                record.seq = seq_protein
            seqs.append(record)

if args['translate_table']!=0:
    print('sequences were translated with the NCBI codon table:', args['translate_table'])
    
if len(not_multiple_of_3) > 0:
    print('summary: sequences not multiple of 3:')
    print(df_id.loc[df_id['gene_id'].isin(not_multiple_of_3),:])


# In[5]:


print('# output sequences:', len(seqs))
with open(args['file_out'], "w") as output_handle:
    Bio.SeqIO.write(seqs, output_handle, "fasta")
print('Completed: ensembl_get_longest_transcript')

