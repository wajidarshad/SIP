# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 16:19:12 2022

@author: Dr. Wajid A. Abbasi
"""

from Bio import SeqIO
from SIP import *

    
    
all_seqs_recs=list(SeqIO.parse("example.fasta", "fasta"))#Give your fasta file instead of example.fasta

for rec in all_seqs_recs:
    print (rec.id)
    if validate(rec.seq)==0:
        print('Input sequence is not valid')
    #print ('Interaction Prediction Score:',predict_CaM_interaction(rec.seq))
    print ('Binding Site Prediction:')
    binding_site_results=predict_BS(rec.seq)
    plot_binding_site(binding_site_results, rec.id)