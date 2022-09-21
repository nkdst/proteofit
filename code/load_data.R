
## libraries

## internal variables
.READCOUNT.GASTROC.PATH <- "./data/readcount_genename_gastroc.txt"
.READCOUNT.SOLEUS.PATH <- "./data/readcount_genename_soleus.txt"

counts.gastroc <- read.delim(.READCOUNT.GASTROC.PATH)
counts.soleus <- read.delim(.READCOUNT.SOLEUS.PATH)
