#! /usr/bin/env python
import sys, os;
import numpy as np;
from cvxopt import matrix;
from cvxopt import solvers;

def cvxopt_run(file = 'test.npz'):

    params = np.load(file);
    P = matrix(params['p']);
    q = matrix(pearams['q']);    
    G = matrix(np.array([[-1,0],[0,-1],[-1,-3],[2,5],[3,4]]), tc='d');
    h = matrix(np.array([0,0,-15,100,80]), tc='d');
    sol = solvers.qp(P, q, G, h);

    print sol['x'];
    print sol['primal objective'];


def main(argv):
    cvxopt_run();
    return;

if __name__ == "__main__":
    main(sys.argv);
