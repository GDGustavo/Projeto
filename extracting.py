from __future__ import division
import matplotlib
import warnings
import numpy as np
from pylab import * 
from numpy import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.text as text
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm
from matplotlib.colors import LightSource
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.offsetbox import AnchoredOffsetbox
from matplotlib.patches import FancyArrowPatch

a2=loadtxt('Energy1.txt')
a2 = a2 + a2[1] - a2[0] - a2[9];
print(a2)
