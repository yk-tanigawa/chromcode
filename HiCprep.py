#! /usr/bin/env python
import sys, os;

class Hic():
    def __init__(self, dirname, chr, res):
        def chr_str2num(str):
            return int(str[3:]);
        def res_str2num(str):
            if(str == '1kb'):
                return 1000;
            elif(str == '5kb'):
                return 5000;
            elif(str == '10kb'):
                return 10000;
            elif(str == '25kb'):
                return 25000;
            elif(str == '50kb'):
                return 50000;
            elif(str == '100kb'):
                return 100000;
            elif(str == '250kb'):
                return 250000;
            elif(str == '500kb'):
                return 500000;
            elif(str == '1Mb'):
                return 1000000;
            else:
                return None;
        self.dirname = dirname;
        self.chrstr = chr;
        self.resstr = res;
        self.chrnum = chr_str2num(chr);
        self.resnum = res_str2num(res);

        # load normalize vector and expected value vector
        self.KRnorm = [];
        self.VCnorm = [];
        self.SQRTVCnorm = [];
        self.RAWexpected = [];
        self.KRexpected = [];
        self.VCexpected = [];
        self.SQRTVCexpected = [];
        with open(dirname + '/' + chr + '_' + res + '.KRnorm', 
                  'r') as fKRnorm:
            for l in fKRnorm:
                self.KRnorm.append(float(l[:-1]));
        with open(dirname + '/' + chr + '_' + res + '.VCnorm', 
                  'r') as fVCnorm:
            for l in fVCnorm:
                self.VCnorm.append(float(l[:-1]));
        with open(dirname + '/' + chr + '_' + res + '.SQRTVCnorm', 
                  'r') as fSQRTVCnorm:
            for l in fSQRTVCnorm:
                self.SQRTVCnorm.append(float(l[:-1]));
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
        fname_head = self.dirname + '/' + self.chrstr + '_' + self.resstr;
        if(norm == 'KR'):
            norm_v = self.KRnorm;
            norm_str = norm;
        elif(norm == 'VC'):
            norm_v = self.VCnorm;
            norm_str = norm;
        elif(norm == 'SQRTVC'):
            norm_v = self.SQRTVCnorm;
            norm_str = norm;
        else:
            norm_str = 'noNorm';

        if(exp == 'RAW'):
            exp_v = self.RAWexpected;
            exp_str = exp;
        elif(exp == 'KR'):
            exp_v = self.KRexpected;
            exp_str = exp;
        elif(exp == 'VC'):
            exp_v = self.VCexpected;
            exp_str = exp;
        elif(exp == 'SQRTVC'):
            exp_v = self.SQRTVCexpected;
            exp_str = exp;
        else:
            exp_str = 'noOE'

        with open(fname_head + '.RAWobserved', 
                  'r') as f:
            wfname = '.'.join([fname_head, norm_str, exp_str, 'dat']);
            with open(wfname, 'w') as out:
                [i, j, mij] = [0, 0, 0.0];
                for l in f:
                    (si, sj, smij) = l[:-1].split();
                    [i, j, mij] = [int(si), int(sj), float(smij)];
                    if(((min == None) or (min <= abs(i - j))) and
                       ((max == None) or (abs(i - j) <= max))):
                        if (norm != None):
                            mij /= (norm_v[i / self.resnum] * norm_v[j / self.resnum]);
                        if (exp != None):
                            mij /= exp_v[abs(i - j) / self.resnum];
                        out.write('{0} {1} {2}'.format(i, j, mij));
                        #print '{0} {1} {2}'.format(i, j, mij);


def hic_prep(datadir, chr, res):
    hic.extract(max = 1000000, norm = None, exp = 'RAW');
    hic.extract(max = 1000000, norm = 'KR', exp = 'KR');
    hic.extract(max = 1000000, norm = 'VC', exp = 'VC');
    hic.extract(max = 1000000, norm = 'SQRTVC', exp = 'SQRTVC');

def main(argv):    
    datadir='../../data/GM12878_conbined_1kb_intra_chr21_MAPQGE30';
    chr = 'chr21';
    res = '1kb';
    hic_prep(datadir, chr, res);
    return;

if __name__ == "__main__":
    main(sys.argv);
