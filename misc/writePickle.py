import job
import sys
import pickle

jobs = []

pickle.dump(jobs,open(sys.argv[1],'wb'))