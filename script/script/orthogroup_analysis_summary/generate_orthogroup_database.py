
# coding: utf-8

# In[1]:


import numpy, pandas, sys, os, re, time, sqlite3, sqlalchemy, datetime, argparse, linecache
from sqlalchemy.types import *

pandas.options.display.max_rows=100
pandas.options.display.max_columns=1000


# In[2]:


# Parameters
# mode = "debug" # "debug" or "batch"
print("Setting parameters...")
if sys.argv[1] == "-f":
    mode = "debug"
else:
    mode = "batch"

if mode=="debug":
    ensembl_dir='/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/'
    params = dict()
    params['overwrite'] = 1
    params['dbpath'] = ensembl_dir+'Ensembl.91.orthogroup.db'
    params['dir_stat_tree'] = ensembl_dir+'orthogroup/stat.tree'
    params['dir_stat_branch'] = ensembl_dir+'orthogroup/stat.branch'
elif mode=="batch":
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', metavar='bool', default=0, type=bool, help='')
    parser.add_argument('--dbpath', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--dir_stat_tree', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--dir_stat_branch', metavar='PATH', default='', type=str, help='')    
    args = parser.parse_args()  
    params = dict()
    for attr in [a for a in dir(args) if not a.startswith('_')]:
        params[attr] = getattr(args, attr)

dir_db = os.path.dirname(params['dbpath'])+'/'
file_db = os.path.basename(params['dbpath'])


# In[3]:


print(datetime.datetime.today(), ': start : database connection')
start = time.time()
os.chdir(dir_db)
if (params['overwrite']) & (os.path.exists(file_db)):
    os.remove(file_db)

conn = sqlalchemy.create_engine("sqlite:///"+file_db)
print(datetime.datetime.today(), ': end : database connection')


# In[4]:


tables = pandas.read_sql_query(sql="SELECT name FROM sqlite_master WHERE type='table'", con=conn)['name'].values
print('Existing tables:', tables)


# In[5]:


print('Sorting input files according to column numbers.', flush=True)

infiles = dict()
columns = dict()
num_columns = dict()
indirs = {'tree':params['dir_stat_tree'],'branch':params['dir_stat_branch'],}
for stat in indirs.keys():
    infiles[stat] = os.listdir(indirs[stat])
    print(stat, ': number of infiles =', len(infiles[stat]), flush=True)
    num_columns[stat] = list()
    for infile in infiles[stat]:
        file_path = os.path.join(indirs[stat], infile)
        num_columns[stat].append(linecache.getline(file_path, 1).count('\t')+1)
    infiles[stat] = (pandas.Series(infiles[stat])[numpy.argsort(num_columns[stat])[::-1]]).tolist()
    df_maxcol = pandas.read_csv(os.path.join(indirs[stat], infiles[stat][0]), sep="\t", header=0, index_col=False)
    columns[stat] = df_maxcol.columns
    print(stat, ': num_col =', len(columns[stat]), ': columns =', columns[stat], flush=True)


# In[6]:


print(datetime.datetime.today(), ': start : add infile to database', flush=True)

for stat in indirs.keys():
    for i in range(len(infiles[stat])):
        og = re.sub('\..*', '', infiles[stat][i])
        file_path = os.path.join(indirs[stat], infiles[stat][i])
        if (os.path.getsize(file_path)==0):
            print('WARNING: skipped due to file size = 0:', infiles[i], flush=True)
        else:
            df = pandas.read_csv(file_path, sep="\t", header=0, index_col=False, dtype=None, engine="c", na_filter=True, low_memory=True, chunksize=None)
            if (all([ c in columns[stat] for c in df.columns ])):
                df["orthogroup"] = og
                df.to_sql(name=stat, con=conn, if_exists="append", index=False, dtype=None, chunksize=1024*8)
            else:
                no_match_columns = df.columns[~pandas.Series([ c in columns[stat] for c in df.columns ])]
                print('WARNING: Column names do not match to the database. Skipped:', infiles[stat][i], flush=True)
                print(no_match_columns, flush=True)
        if (i%100==0):
            print(datetime.datetime.today(), ':', stat, ': processing', i, 'th file.', flush=True)

print(datetime.datetime.today(), ': end : add infile to database', flush=True)


# In[7]:


tables = pandas.read_sql_query(sql="SELECT name FROM sqlite_master WHERE type='table'", con=conn)['name'].values
print('Existing tables:', tables)
db_treeids = dict()
for table in tables:
    try:
        db_treeids[table] = pandas.read_sql_query(sql="SELECT DISTINCT orthogroup FROM "+table, con=conn, index_col=None, coerce_float=True, chunksize=None)
        db_treeids[table] = set([ item for sublist in db_treeids[table].values for item in sublist ])
        print('# tree_id in the database table -', table, ':', len(db_treeids[table]))
    except:
        print('Failed to retrieve tree_ids from', table, 'in', file_db)
        db_treeids[table] = set()


# In[8]:


tables = pandas.read_sql_query(sql="SELECT name FROM sqlite_master WHERE type='table'", con=conn)['name'].values
index_columns = dict()
for table in tables:
    print('\n', table)
    sql = "PRAGMA TABLE_INFO("+table+")"
    columns = pandas.read_sql_query(sql=sql, con=conn, index_col=None, coerce_float=True, chunksize=None)
    print(columns.loc[:,['name','type']].values)
    index_columns[table] =columns['name'].values


# In[9]:


conn.dispose()
print(datetime.datetime.today(), ': completed!')

