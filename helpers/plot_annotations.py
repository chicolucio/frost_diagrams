import numpy as np
from matplotlib import pyplot as plt


class Annotations:

    def __init__(self, x_data, y_data, label, color, facecolor, txt_width=None, txt_height=None, xmove=0, ymove=-0.05):
        self.x_data = x_data
        self.y_data = y_data
        self.label = label
        self.color = color
        self.facecolor = facecolor
        self.xmove = xmove
        self.ymove = ymove
        self.ax = plt.gca()

        if (txt_width is None) and (txt_height is None):
            self.txt_height = 0.04 * (self.ax.get_ylim()[1] - self.ax.get_ylim()[0])
            self.txt_width = 0.03 * (self.ax.get_ylim()[1] - self.ax.get_ylim()[0])
        else:
            self.txt_height = txt_height
            self.txt_width = txt_width

    def _get_text_positions(self):
        # adapted from https://stackoverflow.com/questions/8850142/matplotlib-overlapping-annotations/10739207
        a = list(zip(self.y_data, self.x_data))
        text_positions = self.y_data.copy()
        for index, (y, x) in enumerate(a):
            local_text_positions = [i for i in a if i[0] > (y - self.txt_height)
                                    and (abs(i[1] - x) < self.txt_width * 2) and i != (y, x)]
            if local_text_positions:
                sorted_ltp = sorted(local_text_positions)
                if abs(sorted_ltp[0][0] - y) < self.txt_height:  # True == collision
                    differ = np.diff(sorted_ltp, axis=0)
                    a[index] = (sorted_ltp[-1][0] + self.txt_height, a[index][1])
                    text_positions[index] = sorted_ltp[-1][0] + self.txt_height
                    for k, (j, m) in enumerate(differ):
                        # j is the vertical distance between words
                        if j > self.txt_height * 2:  # if True then room to fit a word in
                            a[index] = (sorted_ltp[k][0] + self.txt_height, a[index][1])
                            text_positions[index] = sorted_ltp[k][0] + self.txt_height
                            break
        return text_positions

    def text_plotter(self):
        # adapted from https://stackoverflow.com/questions/8850142/matplotlib-overlapping-annotations/10739207

        for x, y, l, t in zip(self.x_data, self.y_data, self.label, self._get_text_positions()):
            self.ax.text(x + self.txt_width + self.xmove, t + self.ymove, fr'$\mathrm{{{self.label[x]}}}$',
                         rotation=0, color=self.color, fontsize=20,
                         bbox=dict(facecolor=self.facecolor, alpha=0.1))
            if y != t:
                self.ax.arrow(x, t, 0, y - t, color='red', alpha=0.3, width=self.txt_width * 0.1,
                              head_width=self.txt_width, head_length=self.txt_height * 0.5,
                              zorder=0, length_includes_head=True)