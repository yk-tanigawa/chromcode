#! /usr/bin/env python
import sys, os;
import numpy as np;
import matplotlib.pyplot as plt;

# try to develop a visualize function
def draw_heatmap(data, row_labels, column_labels):
    fig, ax = plt.subplots();
    heatmap = ax.pcolor(data, cmap=plt.cm.Blues);

    ax.set_xticks(np.arange(data.shape[0]) + 0.5, minor=False);
    ax.set_yticks(np.arange(data.shape[1]) + 0.5, minor=False);

    ax.invert_yaxis();
    ax.xaxis.tick_top();

    ax.set_xticklabels(row_labels, minor=False);
    ax.set_yticklabels(column_labels, minor=False);
    plt.show();
    # plt.savefig('image.png');

    return heatmap;

def main(argv):
    k = 5;
    a = np.ones((1 << (2 * k), 1 << (2 * k)), dtype = np.float64);
    for i in xrange(1 << (2 * k)):
        for j in xrange(1 << (2 * k)):
            a[i][j] = i * (1 << (2 * k)) + j;
    #print a;
    lab = [''] * (1 << (2 * k));
    labnum = 8;
    for i in xrange(labnum):
        lab[i * (1 << (2 * k)) / labnum] = i * (1 << (2 * k)) / labnum;
    lab[-1] = (1 << (2 * k)) - 1;
    draw_heatmap(a, lab, lab);
    return;

if __name__ == "__main__":
    main(sys.argv);
