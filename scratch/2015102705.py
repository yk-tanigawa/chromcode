#! /usr/bin/env python
import sys, os;

# want to test set to list

def main(argv):
    s = set([1, 2, 3]);
    print s;
    l = list(s);
    l.sort();
    print l;
    return;

if __name__ == "__main__":
    main(sys.argv);
