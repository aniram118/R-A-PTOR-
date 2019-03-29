### Code has been borrowed and modified from lrplots(https://github.com/ahcm/longread_plots)

import matplotlib

matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import datetime
import seaborn as sns
import pandas as pd
import numpy as np
import collections

matplotlib.rcParams['axes.titlesize'] = 13
matplotlib.rcParams['axes.labelsize'] = 12
matplotlib.rcParams['xtick.labelsize'] = 11
matplotlib.rcParams['ytick.labelsize'] = 11
matplotlib.rcParams['figure.titlesize'] = 20

def fastq_reader(filename):
    i = 0
    with open(filename, 'r') as infile:
        name = infile.readline().rstrip()
        while True:
            i += 1
            seq = ""
            for s in infile:
                if s[0] == '+':
                    commentp = s.rstrip()
                    break
                else:
                    seq += s.rstrip()
            qual = ""
            for q in infile:
                if len(qual) > 0 and q[0] == '@':
                    yield name, seq, qual
                    name = q.rstrip()
                    break
                else:
                    qual += q.rstrip()
            else:
                yield name, seq, qual
                return


def get_fastq_info(filename):
    seq_lengths = []
    mean_qualities = []
    kmers_start = []
    kmers_end = []
    nts_A = []
    nts_G = []
    nts_T = []
    nts_C = []
    nts_U = []

    for head, seq, qual in fastq_reader(filename):

        seq_lengths.append(len(seq))
        seq=seq.replace('T','U')
        mean_qualities.append(round(np.mean(bytearray(qual, "ascii")) - 33, 2))

        kmers_start.append(seq[0:4])
        kmers_end.append(seq[-4:])

        ntc = collections.Counter()
        ntc.update(seq)
        nts_A.append(ntc['A'])
        nts_G.append(ntc['G'])
        nts_T.append(ntc['T'])
        nts_C.append(ntc['C'])
        nts_U.append(ntc['U'])

    df = pd.DataFrame({"seq_length": seq_lengths,
                       "mean_quality": mean_qualities,
                       "kmers_start": kmers_start,
                       "kmers_end": kmers_end,
                       "nt_A": nts_A,
                       "nt_G": nts_G,
                       "nt_T": nts_T,
                       "nt_C": nts_C,
                       "nt_U": nts_U,
                       })
    return df

def xlabel_pos_right(ax):
    label = ax.xaxis.get_label()
    x_lab_pos, y_lab_pos = label.get_position()
    label.set_position([1.0, y_lab_pos])
    label.set_horizontalalignment('right')
    ax.xaxis.set_label(label)

## Plots

def plot_nt_content(ax, df):
    ax.set_title("Nucleotide Content", loc="left")
    ax.set_xlabel("NT")
    ax.set_ylabel("Percent")
    nt = {}
    nt['A'] = np.sum(df.nt_A)
    nt['T'] = np.sum(df.nt_T)
    nt['G'] = np.sum(df.nt_G)
    nt['C'] = np.sum(df.nt_C)
    nt['U'] = np.sum(df.nt_U)
    s = np.sum((nt['A'], nt['T'], nt['G'], nt['C'], nt['U']))
    for k, v in nt.items(): nt[k] = (v / s) * 100
    if nt['U'] > 0:
        sns.barplot(['G', 'C', 'A', 'U'], [nt['G'], nt['C'], nt['A'], nt['U']], ax=ax, palette='dark')
    else:
        sns.barplot(['G', 'C', 'A', 'T'], [nt['G'], nt['C'], nt['A'], nt['T']], ax=ax, palette='dark')

def plot_length_binned(ax, df, num_bins=1000):
    ax.set_title("Histogram of sequence lengths 1000 bins", loc="left")
    ax.set_xlabel("Lengths")
    ax.set_ylabel("Count (Log)")
    ax.set_xlim(left=0, right=6000)
    ax.hist(df.seq_length, num_bins,color='green')
    ax.set_yscale('log', nonposy='clip')
    xlabel_pos_right(ax)

def plot_qualities(ax, df):
    ax.set_title("Histogram of Mean Qualities", loc="left")
    ax.hist(df.mean_quality, 160, range=(0, 40), color='green')  # , s=point_size) #, range=(0,20000))
    ax.set_xlim(left=0, right=40)
    ax.set_xlabel("Quality")
    xlabel_pos_right(ax)

def plot_hist2d_helper(ax, xbins, ybins, x, y):
    counts, _, _ = np.histogram2d(x, y, bins=(xbins, ybins))
    ax.pcolormesh(xbins, ybins, counts.T, cmap="Greens", norm=matplotlib.colors.LogNorm())

def plot_quality_seqlength(ax, df):
    ax.set_title("Read Length vs Quality", loc="left")
    max_len = np.max(df.seq_length)
    xbins = np.arange(0, 40, 0.1)
    ybins = np.arange(1, np.log(max_len), 0.01) ** 10
    plot_hist2d_helper(ax, xbins, ybins, df.mean_quality, df.seq_length)
    ax.set_xlim(left=0, right=40)
    ax.set_xlabel("Mean Quality")
    ax.set_ylabel("Read Length Log10")
    xlabel_pos_right(ax)
    ax.set_yscale('log', nonposy='clip')

# Assemble Plots
def plot_fastq_info(df,filename):
    (fig, axsa) = plt.subplots(ncols=2, nrows=2)
    axs = axsa.flatten().tolist()
    axs.reverse()
    point_size = 8
    marker_size = 3

    plot_length_binned(axs.pop(), df, 1000)

    plot_qualities(axs.pop(), df)

    plot_quality_seqlength(axs.pop(), df)

    plot_nt_content(axs.pop(), df)

    plt.subplots_adjust(left=0.2, wspace=1.0, top=2.8)

    plt.suptitle(filename[:-12] + ' ' + '3\' UMR Analysis Plots')
    fig.set_size_inches(21, 12)
    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    return plt

def plotter(filename):
    df_filename = filename
    # print(filename)
    df = get_fastq_info(filename)
    plt = plot_fastq_info(df,filename)
    plt.savefig(filename[:-6] + ".pdf") 
    return plt