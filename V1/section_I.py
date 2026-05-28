import numpy as np
import scipy.integrate 
from scipy.integrate import quad
from scipy.optimize import root, fsolve, least_squares, root_scalar
from scipy.spatial import Delaunay
from matplotlib.path import Path
from xlwings import func

from beton import *   # ← fonctionne en local et sur GitHub
from acier import *   # ← fonctionne en local et sur GitHub	

"""
========================================================================
CALCUL DE SECTIONS EN I — BÉTON ARMÉ (EC2)
xlwings Lite 
========================================================================
Conventions :
  - Déformations en ‰ (mm/m)
  - Contraintes en MPa
  - Aires en mm²  (entrée cm² × 100, sortie cm²)
  - Es = 200 000 MPa
  - Axe y = 0 à la fibre inférieure, positif vers le haut
  - Compression positive en béton EC2
========================================================================
"""
