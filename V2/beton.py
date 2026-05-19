# beton.py
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
	return np.where(x <0, -1 ,1)
@func		
def eps_c2(fck_array):
    """Renvoie la déformation eps_c2 selon formule 8.4, en mm/m"""
    # arrondi = faux pour avoir la valeur selon la formule
  
    # Conversion en array si ce n'est pas déjà le cas
    fck = np.asarray(fck_array)
    
    # Application de la loi : 2 si fck <= 50, sinon formule complexe
    # np.where(condition, valeur_si_vrai, valeur_si_faux)
    eps = 2.0
    return eps	
@func
def eps_cu2(fck):
    """Renvoie la déformation eps_cu2 selon formule 8.4 en mm/m"""

    return 3.5
@func    
def eps_n(fck_array):
    """Renvoie le coefficient n selon formule 8.4"""
    # arrondi = faux pour avoir la valeur selon la formule
  
    fck = np.asarray(fck_array)
    
    
    eps = 2
    return eps	
@func		
def sigma_c_n(eps_c_array, n_array):
    """
    Version vectorisée pour la contrainte béton fissuré.
    Module Acier = 200 000 MPa.
    n = coefficient d'équivalence (ex: 15).
    eps_c_array : déformation en mm/m.
    """
    # On s'assure que ce sont des tableaux NumPy
    eps = np.asarray(eps_c_array)
    n = np.asarray(n_array)

    # sigma = (Es / n) * eps
    # Le /1000 convertit les mm/m en déformation relative (sans unité)
    sigma = (200000 / n) * (eps / 1000)

    # Le béton fissuré ne reprend pas de traction (eps < 0)
    # np.maximum(0, sigma) remplace toutes les valeurs négatives par 0
    return np.maximum(0, sigma)
@func
def eps_c_n(sig_c_array, n_array):
    """
    Version vectorisée : renvoie la déformation eps_c (en mm/m)
    sig_c_array : contrainte béton (MPa)
    n_array : coefficient d'équivalence Es/Ec
    """
    sig = np.asarray(sig_c_array)
    n = np.asarray(n_array)

    # Calcul de la déformation : epsilon = (n / Es) * sigma * 1000
    # Le * 1000 permet de repasser en mm/m (unité standard Eurocode)
    epsilon = (1000 * n / 200000) * sig

    # Gestion du cas sig_c <= 0 (Béton fissuré/tendu)
    # On remplace les valeurs où sig <= 0 par np.nan (ou 0.0 selon votre besoin)
    return np.where(sig > 0, epsilon, 0.0)
@func		
def sigma_c_pararect(fck, fcd, eps_c_array):
    """
    Version vectorisée de la loi parabole-rectangle (EC2, 3.1.7).
    fck, fcd : valeurs scalaires (ou tableaux)
    eps_c_array : tableau des déformations (mm/m)
    """
    eps_c = np.asarray(eps_c_array)
    
    # Récupération des paramètres (en utilisant nos versions vectorisées précédentes)
    e2 = eps_c2(fck)
    n = eps_n(fck)

    # Définition des conditions
    conditions = [
        (eps_c <= 0),                          # 1. Zone tendue (fissurée)
        (eps_c > 0) & (eps_c <= e2),           # 2. Zone parabolique
        (eps_c > e2)                           # 3. Zone rectangulaire (palier)
    ]

    # Définition des formules correspondantes
    choix = [
        0.0,                                   # 1. sigma = 0
        fcd * (1 - (1 - eps_c / e2)**n),       # 2. sigma = fcd * [1 - (1 - eps/eps2)^n]
        fcd                                    # 3. sigma = fcd
    ]

    return np.select(conditions, choix)
@func
def epsilon_c_pararect(fck, fcd, sig_c_array):
    sig = np.asarray(sig_c_array)
    e2 = eps_c2(fck)
    n = eps_n(fck)

    # On sature pour éviter les erreurs de racine (1/n)
    sig_safe = np.clip(sig, 0, fcd * 0.99999999999999999999999)

    # Calcul analytique
    epsilon = e2 * (1 - (1 - sig_safe / fcd)**(1/n))

    # Gestion des masques
    # Si sigma <= 0 : déformation = 0.0
    # Si sigma >= fcd : déformation = e2 (début du palier)
    resultat = np.select([sig <= 0, sig >= fcd], [0.0, e2], default=epsilon)
    
    return resultat
@func
def sigma_c_n1(eps_c, n):
    """
    Loi béton fissuré. 
    Utilise une fonction de transfert pour annuler la traction sans np.where.
    """
    eps = np.asarray(eps_c)
    eps_abs = np.abs(eps)
    
    # Le terme (eps + eps_abs) / (2 * eps_abs + 1e-15) agit comme un filtre :
    # Si eps > 0  => (eps + eps) / (2 * eps) = 1
    # Si eps < 0  => (eps - eps) / (2 * eps) = 0
    # Le 1e-15 est le "garde-fou" contre la division par zéro.
    
    filtre = (eps + eps_abs) / (2 * eps_abs + 1e-30)
    sigma = (200000 / n * eps / 1000) * filtre
    return sigma



@func
def sigma_c_pararect1(fck, fcd, eps_c_array):
    eps_c=np.array(eps_c_array)
    eps_c2_val = eps_c2(fck)
    eps_n_val = eps_n(fck)
    term1 = np.abs(1 - eps_c / eps_c2_val)
    term2 = (1 - eps_c / eps_c2_val)
    sigma = fcd * (1 - (term1 + term2) / (2 * term1 + 1e-30) * term1 ** eps_n_val) * (eps_c + np.abs(eps_c)) / (2 * np.abs(eps_c) + 1e-30)
    return sigma
