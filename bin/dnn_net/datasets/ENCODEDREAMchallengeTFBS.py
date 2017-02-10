from __future__ import absolute_import
import os
import cPickle
from glob import glob
import gzip

import random
import numpy as np
from keras.utils import np_utils

from ..utils.lambdas import myself, flatten, nuc2ord, listComplement, label2int
from ..utils.data_utils import get_meta, get_label_count, get_label, get_bam5p, get_seq

""" dta_dir should contain subfolders called misc, labels, & DNAse"""
"""N is set to the Nanog's B label count"""
def load_data(tf, N=65836, train_frac=0.8, version = 0.1, remove_cell=['SK-N-SH'], \
              dta_dir='/home/biter/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484', \
              force_read=False, verbose=1):

    pkl_fle = 'v' + str(version) + '.' + tf + '.data.pkl.gz'
    try:
        if force_read:
            raise Exception('forced_read')
        with gzip.open(pkl_fle, 'rb') as handle:
            (X_train, y_train), (X_test, y_test) = cPickle.load(handle)
        if verbose > 1: print "Read " + pkl_fle
    except:
        meta_fle   = os.path.join(dta_dir, 'metadata')
        genome_fle = os.path.join(dta_dir, 'misc', 'hg19.genome.fa.gz')
        label_fle  = glob(os.path.join(dta_dir, 'labels', '*' + tf + '*train.labels.tsv.gz'))[0]
        DNAse_dir  = glob(os.path.join(dta_dir, 'DNAse'))[0]

        count_dic, Nx, Nnew, genomic_window_size = \
            get_label_count(label_fle, N, remove_cell=remove_cell, force_read=force_read, verbose=verbose)

        label_dic = get_label(label_fle, count_dic, genomic_window_size, Nx, \
                              force_read=force_read, verbose=verbose)
        bam5p_dic = get_bam5p(DNAse_dir, label_dic, fle_tag=tf, force_read=force_read, verbose=verbose)
        seq_dic   = get_seq(genome_fle, label_dic, fle_tag=tf, force_read=force_read, verbose=verbose)
    #    shape_dic = get_shape(shape_fle, label_dic, fle_tag=tf, force_read=force_read, verbose=verbose)

        features = 5 # TODO set
        Nnew2 = Nnew // 2
        X = {'B':np.zeros((Nnew2, genomic_window_size * features), dtype=np.float),
             'U':np.zeros((Nnew2, genomic_window_size * features), dtype=np.float)}
        y = {'B':np.zeros((Nnew2), dtype=np.int),
             'U':np.zeros((Nnew2), dtype=np.int)}
        i = {'B':-1, 'U':-1}
        # save space by deleting
        for cell in seq_dic.keys():
            for label in seq_dic[cell].keys():
                for chrom in seq_dic[cell][label].keys():
                    for ni in range(len(bam5p_dic[cell][label][chrom])):
                        i[label] += 1
                        y[label][i[label]] = label2int(label)
                        #if cell == 'K562' and chrom == 10: 
                        #    print label, label2int(label)
                        for gi in range(genomic_window_size):
                            i1 = gi * features
                            i2 = i1 + features
                            nuc = seq_dic[cell][label][chrom][ni][gi]
                            if not nuc in 'ACGT':
                                nuc = random.choice(['A','C','G','T'])
                            X[label][i[label]][i1:i2] = [ bam5p_dic[cell][label][chrom][ni][gi] ] + nuc2ord(nuc)            
                    del label_dic[cell][label][chrom]
                    del bam5p_dic[cell][label][chrom]
                    del seq_dic[cell][label][chrom]
                if (verbose > 1):
                    print cell, label, i[label], X[label][i[label]][0:10]
                del label_dic[cell][label]
                del bam5p_dic[cell][label]
                del seq_dic[cell][label]
            del label_dic[cell]
            del bam5p_dic[cell]
            del seq_dic[cell]

        indices = np.arange(Nnew2)
        train_indices = np.sort(np.random.choice(indices, int(train_frac * Nnew2), replace=False))
        test_indices  = listComplement(indices, train_indices)

        X_train = np.append(X['B'][train_indices,:], X['U'][train_indices,:], axis = 0)
        y_train = np.append(y['B'][train_indices], y['U'][train_indices])

        X_test = np.append(X['B'][test_indices,:], X['U'][test_indices,:], axis = 0)
        y_test = np.append(y['B'][test_indices], y['U'][test_indices])

        if verbose > 1: print "Parsed files for " + tf
        with gzip.open(pkl_fle, 'wb') as handle:
            cPickle.dump([(X_train, y_train), (X_test, y_test)], handle, -1)  # -1 is for HIGHEST_PROTOCOL
        if verbose > 1: print "Wrote " + pkl_fle

    if verbose > 0:
        print "FUNCTION " + myself() + " DTA:"
        print 'X_train' + str(X_train.shape)
        print 'y_train' + str(y_train.shape)
        print 'X_test' + str(X_test.shape)
        print 'y_test' + str(y_test.shape)

    return (X_train, y_train), (X_test, y_test)

def test_load_data():
    load_data('test', force_read=False, verbose=2)
    #load_data('JUND', force_read=False, verbose=2)

