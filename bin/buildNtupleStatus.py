#! /usr/bin/env python
import sys, pickle 


try:
    with open("processingStatus.pkl",'rb'): production = pickle.load(open("processingStatus.pkl",'rb'))
except IOError: production = {}
try: production[sys.argv[1]].append(sys.argv[2])
except KeyError: production[sys.argv[1]] = [sys.argv[2]] 
pickle.dump(production,open("processingStatus.pkl",'wb'))

