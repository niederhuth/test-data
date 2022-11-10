import os
import sys

functionsfile=re.sub('data.*','scripts/sequence/py/functions.py',os.getcwd())
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions as functions

#GC content
print('Getting GC content for transcripts')
functions.gc123(sys.argv[1],output=sys.argv[2])
