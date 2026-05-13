# Module pour les propriétés et calculs de l'acier
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
@func
def sgn(x):
	# renvoie signe 1 / -1
	return np.where(x <0, -10 ,10)

@func
def sigma_s_palier(fyd, k, eps_uk,eps_ud, eps_s_array):
    """renvoie la contrainte acier (fe500 classe B) selon 3.2.7 de NF EN 1992-1-1"""
    eps_s =np.array(eps_s_array)
    eps_sd = fyd / 200000 * 1000 

    
    sigma1 = 200000 * eps_s / 1000 # if abs(eps_s) <= eps_sd:
   
    sigma2 = sgn(eps_s)* (fyd + (k * fyd - fyd)/(eps_uk - eps_sd) * (np.abs(eps_s) - eps_sd)) #  else:   
    return  np.where(np.abs(eps_s) > eps_sd, sigma2, sigma1)	
	
@func
def sigma_s_palier1(fyd, k, eps_uk,eps_ud, eps_s_array,a_com):
    """renvoie la contrainte acier (fe500 classe B) selon 3.2.7 de NF EN 1992-1-1"""
    eps_s =np.array(eps_s_array)

    
    sigma1 =sigma_s_palier(fyd, k, eps_uk,eps_ud, eps_s)*a_com # if eps_s > 0.0:
    sigma2 = sigma_s_palier(fyd, k, eps_uk,eps_ud, eps_s)   # else:
    return  np.where(eps_s > 0, sigma1, sigma2)
@func
def eps_s_palier(fyd, k, eps_uk, eps_ud, sig_s_array):
    """
    Version vectorisée analytique : calcule la déformation eps_s à partir de sig_s.
    sig_s_array : tableau des contraintes (MPa)
    """
    sig = np.asarray(sig_s_array)
    Es = 200000
    eps_sd = (fyd / Es) * 1000
    
    # Préparation des paramètres de la pente plastique
    # Attention : si k=1 (palier horizontal), la déformation n'est pas unique.
    # On gère le cas k > 1
    pente_inv = (eps_uk - eps_sd) / (k * fyd - fyd) if k > 1 else 0
    
    abs_sig = np.abs(sig)
    
    # 1. Calcul pour la zone élastique
    eps_elaste = (sig / Es) * 1000
    
    # 2. Calcul pour la zone plastique (branche inclinée)
    eps_plast = np.sign(sig) * (eps_sd + (abs_sig - fyd) * pente_inv)
    
    # 3. Assemblage des conditions
    # Note : Si sig_s == 0, on renvoie 0.0 au lieu de None pour la compatibilité
    conditions = [
        (sig == 0),                          # Cas sigma nul
        (abs_sig <= fyd),                    # Zone élastique
        (abs_sig > fyd) & (abs_sig <= k*fyd) # Zone plastique
    ]
    
    choix = [
        0.0,
        eps_elaste,
        eps_plast
    ]
    
    # On renvoie 0.0 par défaut si sigma est hors limites
    return np.select(conditions, choix, default=0.0)
	
@func
def sigma_s_lin(eps_s_array):
    """
    Version vectorisée de la loi acier linéaire infinie.
    E = 200 000 MPa (pas de palier plastique).
    eps_s_array : déformation en mm/m.
    """
    # On s'assure que l'entrée est un tableau
    eps = np.asarray(eps_s_array)
    
    # sigma = E * epsilon
    # Le /1000 convertit les mm/m en déformation relative (sans unité)
    return 200000 * eps / 1000
	
@func
def eps_s_lin(sig_s_array):
    """
    Version vectorisée : renvoie la déformation eps_s (en mm/m).
    Loi élastique linéaire infinie (E = 200 000 MPa).
    sig_s_array : contrainte acier (MPa).
    """
    # Conversion en array NumPy
    sig = np.asarray(sig_s_array)
    
    # epsilon = (sigma / E) * 1000 pour avoir des mm/m
    return 1000 * sig / 200000	
	
@func
def sigma_s_lin1(eps_s, a_com):
    """
    Loi acier linéaire avec coefficient de pondération a_com.
    Appliqué si eps_s > 0 (traction).
    """
    eps = np.asarray(eps_s)
    # Calcul de la contrainte de base (E = 200 000 MPa)
    sig_base = 200000 * eps / 1000
    
    # On multiplie par a_com uniquement là où eps_s > 0
    return np.where(eps > 0.0, sig_base * a_com, sig_base)

