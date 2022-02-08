
# coding: utf-8

# In[1]:


"""
Usage: inframe_Nstop2gap.py [infasta] [outfasta]
"""


# In[2]:


import Bio.SeqIO, Bio.Seq, re, os, sys, shutil


# In[3]:


# Parameters
# mode = "debug" # "debug" or "batch"
print("Setting parameters...")
if sys.argv[1] == "-f":
    mode = "debug"
else:
    mode = "batch"

if mode=="debug":
    og ='OG0005190'
    wd = '/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/orthogroup/tmp/'+og+'/'
    os.chdir(wd)
    infasta = wd+og+".cds.fasta"
    outfasta = wd+og+".cds.Nstop_masked.fasta"
    codon_table = "Standard"
elif mode=="batch":
    infasta = sys.argv[1]
    outfasta = sys.argv[2]
    codon_table  = sys.argv[3]


# In[13]:


aln = Bio.SeqIO.parse(infasta, "fasta")
aln = list(aln)

for record in aln:
    nucseq = str(record.seq)
    nucseq_degap = nucseq.replace('-', '')
    nucseq_degap_len = len(nucseq_degap)
    if (nucseq_degap_len%3!=0):
        print('removing partial codons:', record.name, 'length:', nucseq_degap_len)
        remove_len = len(nucseq_degap)%3
        remove_chars = nucseq_degap[(nucseq_degap_len-remove_len):nucseq_degap_len]
        for rc in remove_chars:
            pos_replace = nucseq.rfind(rc)
            nucseq = nucseq[:(pos_replace)]+'-'+nucseq[(pos_replace+1):]
        record.seq = Bio.Seq.Seq(nucseq)
    aaseq = record.seq.translate(table=codon_table, to_stop=False, gap="-")
    for match in re.finditer("X+", str(aaseq)):
        num_X = match.end() - match.start()
        nucseq = nucseq[:match.start() * 3]+"-" * num_X * 3+nucseq[match.end() * 3:]
    for match in re.finditer("\*+", str(aaseq)):
        num_stop = match.end() - match.start()
        nucseq = nucseq[:match.start() * 3]+"-" * num_stop * 3+nucseq[match.end() * 3:]
    record.seq = Bio.Seq.Seq(nucseq)

Bio.SeqIO.write(aln, outfasta, "fasta")

