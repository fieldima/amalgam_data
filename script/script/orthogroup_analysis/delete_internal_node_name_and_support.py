
# coding: utf-8

# In[4]:

"""
Usage:
delete_internal_node_name_and_support.py in.nwk out.nwk
"""

import sys
from ete3 import TreeNode

t = TreeNode(newick=sys.argv[1], format=0)
t.write(outfile=sys.argv[2], format=5)


# In[ ]:



