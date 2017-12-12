from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division

import numpy as np
import matplotlib.pyplot as plt

dataEnergy = np.loadtxt("../averageBeamEnergy.out", dtype=None)
dataTime = np.loadtxt("../averageBeamTime.out", dtype=None)

timeFitX = dataTime[20:40, 0]
timeFitY = dataTime[20:40, 1]

timeFit = np.polyfit(timeFitX, timeFitY, 0)
print(timeFit)
timeFitFunc = np.poly1d(timeFit)

plt.subplot(211)
plt.plot(dataEnergy[:, 0], dataEnergy[:, 1])
plt.subplot(212)
plt.plot(timeFitX, timeFitY)
plt.plot(timeFitX, timeFitFunc(timeFitX))

# plt.xlim([20, 80])
# plt.ylim([5800, 6200])

plt.show()