#! /usr/bin/env python
import sys, os;
import numpy as np;
from cvxopt import matrix;
from cvxopt import solvers;

#def cvxopt_run(infile = 'cvxopt_params.npz', outfile = 'cvxopt.npz'):
def cvxopt_run(infile, outfile):
    try:
        params = np.load(infile);
        P = matrix(params['P']);
        q = matrix(params['q']);
    except IOError:
        print 'file not found';
        raise;
    except KeyError:
        print 'parameters P and q are required to solve QP';
        raise;
    else:
        try:
            G = matrix(params['G']);
            h = matrix(params['h']);
        except KeyError:
            sol = solvers.qp(P, q);
        else:
            try:
                A = matrix(params['A']);
                b = matrix(params['b']);
            except KeyError:
                sol = solvers.qp(P, q, G, h);
            else:
                sol = solvers.qp(P, q, G, h, A, b);
        finally:
            print sol['x'];
            print sol['primal objective'];
            x = np.array([i for i in sol['x']]);
            print x;
            np.savez(outfile, x = x);

def main(argv):
    try:
        infile = argv[1];
        outfile = argv[2];
    except IndexError:
        print 'usage: {pname} <input file> <output file>'.format(pname = argv[0]);

        infile = 'cvxopt_params.npz';
        outfile = 'cvxopt.npz';
        print 'example: {pname} {inf} {outf}'.format(pname = argv[0],
                                                   inf = infile,
                                                   outf = outfile);
    else:
        cvxopt_run(infile, outfile);
        return;

if __name__ == "__main__":
    main(sys.argv);
