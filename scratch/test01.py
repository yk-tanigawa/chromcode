#! /usr/bin/env python

import sys;
import subprocess;
import numpy as np;
from collections import deque;
from datetime import datetime;


def test01():
    print subprocess.check_output('ls');

def main(*argv):
    print 1;
    test01();
    return;

if __name__ == "__main__":
    main(sys.argv);
#    if(len(sys.argv) < 4):
#        sys.stderr.write("Some arguments are missing\n");
#        sys.exit(1);
#    else:
#        main(sys.argv);
