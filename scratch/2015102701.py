#! /usr/bin/env python
import sys, os;
import numpy as np;


def numpy_save():
    ndarr1 = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]);
    ndarr2 = np.array([10, 20, 30]);
    np.savez('test.npz', x=ndarr1, y=ndarr2);
    return;

def numpy_load():
    savedata = np.load('test.npz');
    print savedata['x'];
    print savedata['y'];
    


def main(argv):
    numpy_save();
    numpy_load();
    return;

if __name__ == "__main__":
    main(sys.argv);
