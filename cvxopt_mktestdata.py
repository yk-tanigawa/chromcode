#! /usr/bin/env python
import sys, os;
import numpy as np;


def numpy_save(fname):
    P = np.array([[1, 0], [0, 0]], dtype=np.float64);
    q = np.array([3, 4], dtype=np.float64);
    G = np.array([[-1,0],[0,-1],[-1,-3],[2,5],[3,4]], dtype=np.float64);
    h = np.array([0,0,-15,100,80], dtype=np.float64);

    np.savez(fname, P = P, q = q, G = G, h = h);
    return;

def numpy_load():
    savedata = np.load('test.npz');
    print savedata['P'];
    print savedata['q'];
    print savedata['G'];
    print savedata['h'];


def main(argv):
    try:
        fname = argv[1];
    except IndexError:
        print 'usage: {prog} <output file name>'.format(prog = argv[0]);
        fname = 'cvxopt_params.npz';
        print 'usage: {prog} {fname}'.format(prog = argv[0],
                                             fname = fname);
    else:
        numpy_save(fname);
        return;

if __name__ == "__main__":
    main(sys.argv);
