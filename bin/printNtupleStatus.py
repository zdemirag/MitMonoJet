#! /usr/bin/env python
import sys, pickle

production =  pickle.load(open("processingStatus.pkl",'rb'))
for status in production:
    print ''
    print status
    for dataset in production[status]: print dataset
