#! /usr/bin/env python
import sys, os;
import numpy as np;

def pos(hicfile = 'hic.pos.npz', 
               genomefile = 'genome.pos.npz',
               outfile = 'pos.npz'):
    try:
        savedhic = np.load(hicfile);
        savedgen = np.load(genomefile);
        hic_set = set(list(savedhic['hic']));
        gen_set = set(list(savedgen['genome']));
    except IOError:
        print 'file not found';
        raise;
    except KeyError:
        print 'file exists, but parameter not found';
        raise;
    else: 
        np.savez(outfile, 
                 pos = np.array(list(hic_set.intersection(gen_set))));

def main(argv):
    pos();
    return;

if __name__ == "__main__":
    main(sys.argv);
