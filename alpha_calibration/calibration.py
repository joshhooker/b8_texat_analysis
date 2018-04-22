from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import scipy.stats as stats

f = open('calibration.dat', 'w')

data = np.genfromtxt('calibration.in', dtype=None)
y = [3183, 5156, 5485, 5805]

for i in range(len(data)):
  x = [data[i][3], data[i][4], data[i][5], data[i][6]]
  slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
  print data[i][0], data[i][1], data[i][2], slope, intercept
  f.write('%s %s %s %s %s\n' %(data[i][0], data[i][1], data[i][2], slope, intercept))

f.close()