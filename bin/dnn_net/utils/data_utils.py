from __future__ import absolute_import
import cPickle
import glob
import os
import gzip

import numpy as np
from Bio import SeqIO
from pandas import read_table
import pyDNase

from ..utils.lambdas import myself, chrstrip

#------------------FASTA---------------------------
def get_seq(fle, label_dic, fle_tag='TF', force_read=False, verbose=1):

    pkl_fle = fle_tag + '_' + '_'.join(label_dic.keys()) + '.seq.pkl.gz'
    try:
        if force_read:
            raise Exception('forced_read')
        with gzip.open(pkl_fle, 'rb') as handle:
            dta = cPickle.load(handle)
        if verbose > 1: print "Read " + pkl_fle
    except:
        dta = {}
        with gzip.open(fle) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                try:
                    chrom = chrstrip(record.id)
                except:
                    continue #skip non-regular chromosomes
                seq = str(record.seq.upper())
                for cell in label_dic.keys():
                    if not cell in dta.keys():
                        dta[cell] = {}
                    for label in label_dic[cell].keys():
                        if not label in dta[cell].keys():
                            dta[cell][label] = {}
                        if chrom in label_dic[cell][label].keys():
                            Nx = len(label_dic[cell][label][chrom])
                            # strands?
                            #dta[cell][label][chrom] = np.zeros((Nx, window_size * len(strands)), dtype=str)
                            window_size = label_dic[cell][label][chrom][0][1] - label_dic[cell][label][chrom][0][0]
                            # str type???
                            dta[cell][label][chrom] = []
                            gi = 0
                            for grange in label_dic[cell][label][chrom]:
                                i0 = label_dic[cell][label][chrom][gi][0]
                                i1 = label_dic[cell][label][chrom][gi][1]
                                dta[cell][label][chrom].append(seq[i0:i1])
                                gi += 1
                if verbose > 1: print "Parsed chr" + str(chrom)
        if verbose > 1: print "Parsed fasta file " + fle
        with gzip.open(pkl_fle, 'wb') as handle:
            cPickle.dump(dta, handle, -1)  # -1 is for HIGHEST_PROTOCOL
        if verbose > 1: print "Wrote " + pkl_fle

    if verbose > 0:
        print "FUNCTION " + myself() + " DTA:"
        for cell in dta.keys():
            for label in dta[cell].keys():
                i = 0
                printed = False
                for chrom in dta[cell][label].keys():
                    i += len(dta[cell][label][chrom])
                    if not printed:
                        print dta[cell][label][chrom][0][:90]
                        printed = True
                print cell, label, "ALL chroms" + str(i)
        print

    return(dta)

def test_get_seq():
    label_fle = '/home/biter/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484/labels/test.train.labels.tsv.gz'
    label_dic, Nnew = get_label(label_fle, force_read=False, verbose=2)

    #print label_dic
    fle = '/home/biter/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484/misc/hg19.genome.fa.gz'
    get_seq(fle, label_dic, force_read=False, verbose=2)


#------------------BAM-----------------------------
""" returns a dictionary of cellxlabel for average 5p cut"""
""" TODO read filtering in bam file"""
""" TODO: """
""" could not decipher what this function returns; cut counts per position?"""
def get_bam5p(bdir, label_dic, strands=['+', '-'], fle_tag="TF", #genomic_window_size=200,\
              force_read=False, verbose=1):
    pkl_fle = fle_tag + "_" + "_".join(label_dic.keys()) + '.dnase.pkl.gz'
    try:
        if force_read:
            raise Exception('forced_read')
        with gzip.open(pkl_fle, 'rb') as handle:
            dta = cPickle.load(handle)
        if verbose > 1: print "Read " + pkl_fle
    except:
        dta = {}
        for cell in label_dic.keys():
            bam_fles = glob.glob(os.path.join(bdir, '*' + cell + '*.bam'))
            dta[cell] = {}
            for label in label_dic[cell].keys():
                dta[cell][label] = {}
                for chrom in label_dic[cell][label].keys():
                    Nx = len(label_dic[cell][label][chrom])
                    window_size = label_dic[cell][label][chrom][0][1] - label_dic[cell][label][chrom][0][0]
                    #dta[cell][label][chrom] = np.zeros((Nx, window_size * len(strands)), dtype=np.float)
                    dta[cell][label][chrom] = np.zeros((Nx, window_size), dtype=np.float)
                    for bam_fle in bam_fles:
                        reads = pyDNase.BAMHandler(bam_fle, caching=False)
                        gi = 0
                        for grange in label_dic[cell][label][chrom]:
                            # could not decipher what this function returns; cut counts per position?
                            temp = reads["chr%s,%i,%i,+" % (chrom, grange[0], grange[1])]
                            #si = 0
                            for strand in strands:
                                dta[cell][label][chrom][gi,:] += temp[strand]
                                #dta[cell][label][chrom][gi, range(si, si + window_size)] += temp[strand]
                                #si += window_size
                            gi += 1
                    # average per bam_fle for a cell
                    N = 0.0 + len(bam_fles)
                    dta[cell][label][chrom] = dta[cell][label][chrom] / N
            if verbose > 1 : print "Parsed " + str(bam_fles) + ' for average 5p cuts for ' + cell
        with gzip.open(pkl_fle, 'wb') as handle:
            cPickle.dump(dta, handle, -1)  # -1 is for HIGHEST_PROTOCOL
        if verbose > 1: print "Wrote " + pkl_fle

    if verbose > 0:
        print "FUNCTION " + myself() + " DTA:"
        for cell in dta.keys():
            for label in dta[cell].keys():
                i = 0
                for chrom in dta[cell][label].keys():
                    i += len(dta[cell][label][chrom])
                print cell, label, "ALL chroms" + str(i)
        print

    return(dta)

def test_get_bam5p():
    label_fle = '/home/biter/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484/labels/test.train.labels.tsv.gz'
    label_dic, Nnew = get_label(label_fle, force_read=False, verbose=2)

    #print label_dic
    bdir = '/home/biter/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484/DNAse'
    get_bam5p(bdir, label_dic, force_read=False, verbose=2)

#------------------TEXT----------------------------
def get_label_count(fle, N, remove_cell=[], fle_tag='', force_read=False, verbose=1):

    pkl_fle = fle_tag + os.path.basename(fle) + '.count.pkl.gz'
    try:
        if force_read:
            raise Exception('forced_read')
        with gzip.open(pkl_fle, 'rb') as handle:
            count_dic, Nx, Nnew, genomic_window_size = cPickle.load(handle)
        if verbose > 1: print "Read " + pkl_fle
    except:
        with gzip.open(fle, 'rb') as handle:
            header = handle.readline().strip().split('\t')
            count_dic = {}
            for cell in header[3:]:
                if cell in remove_cell:
                    continue
                count_dic[cell] = {}
            for l in handle:
                ls = l.strip().split('\t')
                genomic_window_size = int(ls[2]) - int(ls[1])
                for i in range(3, len(ls)):
                    cell = header[i]
                    if cell in remove_cell:
                        continue
                    if not cell in count_dic.keys():
                        count_dic[cell] = {}
                    label = ls[i]
                    if not label in count_dic[cell].keys():
                        count_dic[cell][label] = 0
                    count_dic[cell][label] += 1

        # set Nx Nnew given label availability
        nb_cells = len(count_dic.keys())
        Nx       = int(N / 2.0 / nb_cells) # pos|neg label per cell are equal
        NxB      = min([ count_dic[i]['B'] for i in count_dic.keys() ])
        NxU      = min([ count_dic[i]['U'] for i in count_dic.keys() ])
        Nx       = min([NxB, NxU, Nx])
        Nnew     = Nx * 2 * nb_cells

        if verbose > 1: print "Parsed " + fle + ' for counts'
        with gzip.open(pkl_fle, 'wb') as handle:
            cPickle.dump([count_dic, Nx, Nnew, genomic_window_size], handle, -1)  # -1 is for HIGHEST_PROTOCOL
        if verbose > 1: print "Wrote " + pkl_fle

    if verbose > 0:
        print "FUNCTION " + myself() + " DTA:"
        print "N=", N, " Nx=", Nx, " Nnew=", Nnew
        print count_dic
        print

    return(count_dic, Nx, Nnew, genomic_window_size)

def test_get_label_count():
    fle = '/home/biter/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484/labels/GABPA.train.labels.tsv.gz'
    fle = '/home/biter/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484/labels/test.train.labels.tsv.gz'

    counts, Nx, Nnew = get_label_count(fle, 100, verbose=2)

def get_label(fle, count_dic, genomic_window_size, Nx, fle_tag='', force_read=False, verbose=1):

    #initially keeps sorted occurance index then granges of sampled label entries
    pkl_fle = fle_tag + os.path.basename(fle) + '.sample.pkl.gz'
    try:
        if force_read:
            raise Exception('forced_read')
        with gzip.open(pkl_fle, 'rb') as handle:
            dta = cPickle.load(handle)
        if verbose > 1: print "Read " + pkl_fle
    except:
        dta = {}
        ite = {} #keeps an iterator on dta
        cnt = {} #keeps current count of labelxcell
        with gzip.open(fle, 'rb') as handle:
            header = handle.readline().strip().split('\t')
            del header[:3]
            for cell in header:
                if not cell in count_dic.keys():
                    continue
                cnt[cell] = dict([('B', 0), ('U', 0)])
                ite[cell] = dict([('B', 0), ('U', 0)])
                dta[cell] = dict([('B', np.zeros((Nx, 3), dtype=np.uint32)),
                                   ('U', np.zeros((Nx, 3), dtype=np.uint32))])
                n = [ count_dic[cell]['B'], count_dic[cell]['U']]
                dta[cell]['B'][:,0] = np.sort(np.random.choice(np.arange(n[0]), Nx, replace=False))
                dta[cell]['U'][:,0] = np.sort(np.random.choice(np.arange(n[1]), Nx, replace=False))
            # select sample granges
            for l in handle:
                ls = l.strip().split('\t')
                grange = [ chrstrip(ls.pop(0)), int(ls.pop(0)), int(ls.pop(0)) ]
                for i in range(len(ls)):
                    cell  = header[i]
                    label = ls[i]
                    if (not cell in count_dic.keys()) or label == 'A':
                        continue
                    itev = ite[cell][label]
                    if itev < Nx and cnt[cell][label] == dta[cell][label][itev, 0]:
                        dta[cell][label][itev] = grange
                        ite[cell][label] = itev + 1
                    cnt[cell][label] += 1

            # put chromosome also as hash key
            for cell in dta.keys():
                for label in dta[cell].keys():
                    temp = dta[cell][label]
                    del dta[cell][label]
                    dta[cell][label] = {}
                    for i in range(temp.shape[0]):
                        chrom = temp[i,0]
                        if chrom in dta[cell][label].keys():
                            dta[cell][label][chrom].append(temp[i,1:])
                        else:
                            dta[cell][label][chrom] = [temp[i,1:]]
        if verbose > 1: print "Parsed " + fle + ' to sample ' + 'granges'

        with gzip.open(pkl_fle, 'wb') as handle:
            cPickle.dump(dta, handle, -1)  # -1 is for HIGHEST_PROTOCOL
        if verbose > 1: print "Wrote " + pkl_fle

    if verbose > 0:
        print "FUNCTION " + myself() + " DTA:"
        for cell in dta.keys():
            for label in dta[cell].keys():
                i = 0
                for chrom in dta[cell][label].keys():
                    i += len(dta[cell][label][chrom])
                print cell, label, "ALL chroms" + str(i)
        print

    return(dta)

def test_get_label():
    fle = '/home/biter/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484/labels/GABPA.train.labels.tsv.gz'
    fle = '/home/biter/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484/labels/test.train.labels.tsv.gz'

    print "will parse original file"
    dta, Nnew, genomic_window_size = get_label(fle, N=9, verbose=2)
    print
    print "will read from pickle file"
    dta, Nnew, genomic_window_size = get_label(fle, N=9, verbose=2)

    print
    print dta

#------------------TEXT----------------------------
""" given a file with tfs as 1st column, returns a dictionary with tf keys """
def get_meta(fle, verbose=1):

    dic = read_table(fle, index_col=0, comment='#').transpose().to_dict()

    if verbose > 0:
        print "FUNCTION " + myself() + " DTA:"
        for i in dic.keys():
            print i, ":", dic[i]
        print

    return(dic)

def test_get_meta():
    meta_file='/home/biter/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484/metadata'
    dic = get_meta(meta_file, verbose=2)

    print dic

