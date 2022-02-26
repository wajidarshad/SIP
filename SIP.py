# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 11:50:26 2022

@author: Dr. Wajid A. Abbasi
"""
from Bio.Data import IUPACData 
import os,random
import numpy as np

window_size=21

def plot_binding_site(BS_results,seq_id):
    if not os.path.exists('binding_site_plots'):
        os.makedirs('binding_site_plots')
    seq_id="".join(x for x in seq_id if x.isalnum())
    import matplotlib.pyplot as plt
    import numpy as np
    BS_plot='binding_site_plots/'+seq_id+'_BS_plot.jpg'
    max_score,max_position,scores=BS_results[0],BS_results[1],BS_results[2]
    fig=plt.figure()
    plt.ioff()# prevent plot to display
    plt.plot(range(1,scores.shape[0]+1),scores,'b-')
    plt.plot(max_position,max_score,'ro')
    plt.xlabel('Position')
    plt.ylabel('Score')
    plt.title('Max score is %f and max position is %d' %(max_score,max_position))
    plt.grid()
    fig.savefig(BS_plot)
def validate(sequence):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWYXBUZ'
    sequence = sequence.upper()
    #if len(set(sequence).difference(amino_acids)):
       # raise ValueError("Input sequence contains non-standard amino acid codes")
    if len(sequence)<21 or len(set(sequence).difference(amino_acids))>0:
        #raise ValueError("Input sequence is shorter than %d amino acids." %window_size)
        return 0
    else:
        return 1

def count_amino_acids(seq): 
         prot_dic = dict((k, 0) for k in IUPACData.protein_letters) 
         for aa in prot_dic: 
              prot_dic[aa] = seq.count(aa) 
         return prot_dic
def one_spectrum_features_PD(seq):#Simply Compute Amino Acid Composition
    features_list=[]
    feature=[]
    AA=["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I","L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    for i in range((len(seq)-window_size)+1):
        wenseq= str(seq[i:window_size+i])
        wenseq=wenseq.replace('U',AA[random.randint(0,len(AA)-1)])
        for j in range(len(wenseq)):
            feature.extend(list(count_amino_acids(wenseq[j]).values()))
        #feature=np.mean(feature,axis=0)
        feature=feature/np.linalg.norm(feature)
        features_list.append(feature)
        feature=[]
    return features_list

def window_scores(features):
    """
    Obtain raw window scores for all windows
    """
    filepath = os.path.join('', 'MIL_1mer-PD_weight_vector.npy')
    W_PD=np.load(filepath)
    scores =  np.dot(features,W_PD)
    max_position = np.argmax(scores)
    max_score = np.round(scores[max_position],decimals=2)
    scores=np.round([np.max(scores[max(0,i-20):min(i+1,len(scores))]) for i in range(len(scores))],decimals=2)
    print ('Position of central residue of maximum scoring window:', max_position+11)
    print ('Maximum Score:', max_score)
    print ('Scores for  central residue of all windows',scores)
    return (max_score,max_position+11,scores)
def predict_BS(sequence):
    if not(validate(sequence)):
        print('Input sequence is shorter than 21 amino acids')
        return
    feats=one_spectrum_features_PD(sequence)
    Scores=window_scores(feats)
    return Scores
