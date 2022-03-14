import numpy as np

from matplotlib.image import NonUniformImage
import matplotlib.pyplot as plt
import pandas as pd
df = pd.read_excel()

hist, bin_edges=np.histogram([1, 2, 1], bins=[0, 1, 2, 3])
hist.sum()

rng = np.random.RandomState(10)  # deterministic random data
a = np.hstack((rng.normal(size=1000), \
               rng.normal(loc=5, scale=2, size=1000)))
_ = plt.hist(a, bins='auto')  # arguments are passed to np.histogram
plt.title("Histogram with 'auto' bins")
plt.show()