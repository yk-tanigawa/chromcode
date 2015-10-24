#! /usr/bin/env python
import sys, os;

class Hic():
    def __init__(self, dirname, chr, res, maxdist = None):
        self.RAWexpected = [];
        self.KRexpected = [];
        self.VCexpected = [];
        self.SQRTVCexpected = [];
        with open(dirname + '/' + chr + '_' + res + '.RAWexpected', 
                  'r') as fRAWexpected:
            for l in fRAWexpected:
                self.RAWexpected.append(float(l[:-1]));
        with open(dirname + '/' + chr + '_' + res + '.KRexpected', 
                  'r') as fKRexpected:
            for l in fKRexpected:
                self.KRexpected.append(float(l[:-1]));
        with open(dirname + '/' + chr + '_' + res + '.VCexpected', 
                  'r') as fVCexpected:
            for l in fVCexpected:
                self.VCexpected.append(float(l[:-1]));
        with open(dirname + '/' + chr + '_' + res + '.SQRTVCexpected', 
                  'r') as fSQRTVCexpected:
            for l in fSQRTVCexpected:
                self.SQRTVCexpected.append(float(l[:-1]));
    def extract(self, min = None, max = None, 
                norm = None, exp = None):
        with open(dirname + '/' + chr + '_' + res + '.RAWobserved', 
                  'r') as f:
            for l in f:
                (si, sj, smij) = l[:-1].split();
                (i, j, mij) = [int(si), int(sj), float(smij)];
                if(abs(i - j) < 1000000):
                    print '{0} {1} {2}'.format(i, j, mij);

    #def __str__(self):
    #    return len(self.RAWexpected);

def openRawdata():
    hicfile='../../data/GM12878_conbined_1kb_intra_chr21_MAPQGE30/chr21_1kb.RAWobserved'
    with open(hicfile, 'r') as f:
        for l in f:
            (si, sj, smij) = l[:-1].split();
            (i, j, mij) = [int(si), int(sj), float(smij)];
            if(abs(i - j) < 1000000):
                print '{0}\t{1}\t{2}'.format(i, j, mij);
    return;
def main(argv):
    datadir='../../data/GM12878_conbined_1kb_intra_chr21_MAPQGE30'
    hic = Hic(datadir, 'chr21', '1kb');
    print len(hic.RAWexpected);
    print len(hic.KRexpected);
    print len(hic.VCexpected);
    print len(hic.SQRTVCexpected);
    
    return;

if __name__ == "__main__":
    main(sys.argv);
