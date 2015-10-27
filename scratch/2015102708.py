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
    a = np.ones((5, 5), dtype = np.float64);
    for i in xrange(5):
        for j in xrange(5):
            a[i][j] = i * 5 + j;
    print a;
    draw_heatmap(a, range(5), range(5));
    return;

if __name__ == "__main__":
    main(sys.argv);
