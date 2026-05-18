# beton.py - V2
#version V2
import datetime as dt
import numpy as np
import pandas as pd
import seaborn as sns
import xlwings as xw
from xlwings import func, script
import matplotlib
import scipy
import math
from scipy.optimize import root, fsolve
from matplotlib.path import Path
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

from beton import *   # ← fonctionne en local et sur GitHub
from acier import *   # ← fonctionne en local et sur GitHub	



# section.py — version finale complète
# Intégration analytique EXACTE par théorème de Green
# ELS et ELU — indépendant du sens de saisie du polygone

