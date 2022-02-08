
# coding: utf-8

# In[1]:


import Bio.SeqIO, Bio.Seq, os, sys


# In[2]:


# Parameters
# mode = "debug" # "debug" or "batch"
print("Setting parameters...")
if sys.argv[1] == "-f":
    mode = "debug"
else:
    mode = "batch"

if mode=="debug":
    wd = '/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/orthogroup/tmp/10_OG0000009_original_0.0001_0.001_0.001/'
    os.chdir(wd)
    infasta = wd+"OG0000009.cds.aln.fasta"
    outfasta = wd+"OG0000009.cds.aln.out.fasta"
elif mode=="batch":
    infasta = sys.argv[1]
    outfasta = sys.argv[2]

# https://stackoverflow.com/questions/32833230/biopython-alignio-valueerror-says-strings-must-be-same-length
# In[3]:


records = Bio.SeqIO.parse(infasta, 'fasta')
records = list(records)
maxlen = max(len(record.seq) for record in records)

for record in records:
    if len(record.seq) != maxlen:
        sequence = str(record.seq).ljust(maxlen, '-')
        record.seq = Bio.Seq.Seq(sequence)
assert all(len(record.seq) == maxlen for record in records)

with open(outfasta, 'w') as f:
    Bio.SeqIO.write(records, f, 'fasta')

