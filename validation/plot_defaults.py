import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib

colors = sns.color_palette(n_colors=6)
sns.set_style('whitegrid')
sns.set_context('paper')
matplotlib.rcParamsDefault['font.size'] = 72
onecolumn_figsize = (8,6)
twocolumn_figsize = (15,6)
fs = 16
legend_fs = fs
tick_fs = 0.8*fs
label_fs = fs
legend_fs = 0.8*fs

markersize = 10
tt_col = colors[0]
lsd_col = colors[1]
beast_col = colors[2]

def shift_point_by_markersize(axes, x, y, markersize):
    """
    Shift overlapping points alonmg x axis by half of the markersize.
    This allows to show better plots with errorbars.
    """
    inv = axes.transData.inverted()
    points = [(i,j) for i,j in zip(x,y)]
    pixels = axes.transData.transform(points)
    res = inv.transform(pixels+(markersize/2,0))
    return res[:,0], res[:, 1]

if __name__ == '__main__':
    pass
