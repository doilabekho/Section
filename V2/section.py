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






@func

def acier_G(p_acier, s_acier):
    """
    Construit la liste des armatures acier [x, y, section],
    en supprimant UNIQUEMENT les lignes non remplies.
    Les sections nulles ou négatives sont conservées.
    """

    if p_acier is None or len(p_acier) == 0:
        return None, None, None

    p_acier = np.atleast_2d(p_acier)
    s_acier = np.atleast_1d(s_acier)

    newacier = []

    nx = min(p_acier.shape[0], s_acier.shape[0])

    for i in range(nx):
        try:
            # nettoyage UNIQUEMENT des coordonnées vides
            if p_acier[i, 0] in (None, "", " ") or p_acier[i, 1] in (None, "", " "):
                continue

            x = float(p_acier[i, 0])
            y = float(p_acier[i, 1])
            s = float(s_acier[i])   # peut être 0 ou négatif → CONSERVÉ ✅

            newacier.append([x, y, s / 10000.0])  # mm² → cm²

        except (ValueError, TypeError, IndexError):
            # ignore uniquement les lignes vraiment invalides
            continue

    return newacier

###########################################################"##"
# travailler avec les liste pour évidement et transformer donnée
#######################################################

@func
def n_liste(data):
    """
    Compte le nombre de séparateurs 'fin / fin',
    insensible à la casse et aux espaces.
    """
    count = 0

    for row in data:
        if len(row) >= 2:
            v0 = str(row[0]).strip().lower()
            v1 = str(row[1]).strip().lower()
            if v0 == "fin" or v1 == "fin":
                count += 1

    return count
@func
def liste_i(data, i):
    """
    Retourne la i-ème sous-liste séparée par 'fin / fin'
    (insensible à la casse et aux espaces).
    """
    fins = []

    for idx, row in enumerate(data):
        if len(row) >= 2:
            v0 = str(row[0]).strip().lower()
            v1 = str(row[1]).strip().lower()
            if v0 == "fin" or v1 == "fin":
                fins.append(idx)

    # sécurité
    if i < 1 or i > len(fins):
        return []

    start = 0 if i == 1 else fins[i - 2] + 1
    end = fins[i - 1]

    return data[start:end]

    return data[start:end]
@func
def points_valides(data):
    pts = []
    for row in data:
        if len(row) >= 2:
            try:
                x = float(row[0])
                y = float(row[1])
                pts.append((x, y))
            except:
                pass
    return pts
@func
def centre_gravite_polygone(data):
    pts = points_valides(data)
    n = len(pts)

    if n < 3:
        return None, None

    A = 0.0
    Cx = 0.0
    Cy = 0.0

    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        cross = x1 * y2 - x2 * y1
        A += cross
        Cx += (x1 + x2) * cross
        Cy += (y1 + y2) * cross

    A *= 0.5
    if A == 0:
        return None, None

    Cx /= (6 * A)
    Cy /= (6 * A)

    return Cx, Cy
@func
def aire_et_centre(data):
    pts = points_valides(data)
    n = len(pts)

    if n < 3:
        return 0.0, None, None

    A = 0.0
    Cx = 0.0
    Cy = 0.0

    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        cross = x1 * y2 - x2 * y1
        A += cross
        Cx += (x1 + x2) * cross
        Cy += (y1 + y2) * cross

    A *= 0.5
    if A == 0:
        return 0.0, None, None

    Cx /= (6 * A)
    Cy /= (6 * A)

    return abs(A), Cx, Cy
@func
def aire_centre_inertie_origine(data):
    pts = points_valides(data)
    n = len(pts)

    if n < 3:
        return 0.0, None, None, 0.0, 0.0

    A = 0.0
    Cx = 0.0
    Cy = 0.0
    Ix = 0.0
    Iy = 0.0

    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        cross = x1 * y2 - x2 * y1

        A += cross
        Cx += (x1 + x2) * cross
        Cy += (y1 + y2) * cross

        Ix += (y1**2 + y1*y2 + y2**2) * cross
        Iy += (x1**2 + x1*x2 + x2**2) * cross

    A *= 0.5
    if A == 0:
        return 0.0, None, None, 0.0, 0.0

    Cx /= (6 * A)
    Cy /= (6 * A)
    Ix /= 12
    Iy /= 12

    return np.abs(A), Cx, Cy, np.abs(Ix), np.abs(Iy)
@func
def caracteristiques_mecaniques(contour, evidement):
    A, Cx, Cy, Ix0, Iy0 = aire_centre_inertie_origine(contour)

    if A == 0:
        return 0.0, None, None, 0.0, 0.0

    A_tot = A
    Cx_tot = A * Cx
    Cy_tot = A * Cy
    Ix_tot = Ix0
    Iy_tot = Iy0

    N = n_liste(evidement)

    for i in range(1, N + 1):
        Li = liste_i(evidement, i)
        Ai, cxi, cyi, ixi, iyi = aire_centre_inertie_origine(Li)

        if Ai > 0:
            A_tot -= Ai
            Cx_tot -= Ai * cxi
            Cy_tot -= Ai * cyi
            Ix_tot -= ixi
            Iy_tot -= iyi

    if A_tot == 0:
        return 0.0, None, None, 0.0, 0.0

    CxG = Cx_tot / A_tot
    CyG = Cy_tot / A_tot

    IxG = Ix_tot - A_tot * CyG**2
    IyG = Iy_tot - A_tot * CxG**2

    return np.abs(A_tot), CxG, CyG, np.abs(IxG), np.abs(IyG)
@func
def translation_points(data, Cx, Cy):
    """
    Transforme une liste de points dans le repère centré en (Cx, Cy)
    """
    pts_trans = []

    for row in data:
        if len(row) >= 2:
            try:
                x = float(row[0])
                y = float(row[1])
                pts_trans.append((x - Cx, y - Cy))
            except:
                pass  # ignore lignes vides ou non numériques

    return pts_trans

@func
def transformation_repere_cg(contour, evidement):
    """
    Calcule le centre de gravité global puis transforme les coordonnées
    du contour et des évidements dans le repère centré (Cx, Cy).
    Retour :
    - contour_cg        : liste de points transformés
    - evidements_cg     : liste de listes (vide si N = 0)
    - Cx, Cy            : centre de gravité global
    """

    # --- Calcul du centre de gravité global
    A, Cx, Cy, Ix, Iy = caracteristiques_mecaniques(contour, evidement)

    if Cx is None or Cy is None:
        return None, [], None, None

    # --- Contour transformé
    contour_cg = translation_points(contour, Cx, Cy)

    # --- Évidements transformés (séparés)
    evidements_cg = []
    N = n_liste(evidement)

    if N == 0:
        # Aucun évidement
        return contour_cg, evidements_cg, Cx, Cy

    for i in range(1, N + 1):
        Li = liste_i(evidement, i)
        Li_cg = translation_points(Li, Cx, Cy)
        evidements_cg.append(Li_cg)

    return  contour_cg, evidements_cg, Cx, Cy

@func
def nettoyer_polygone(points, tol=1e-9):
    """
    Supprime le dernier point s'il est égal au premier
    à une tolérance numérique près.
    """
    if len(points) < 2:
        return points

    x1, y1 = points[0]
    x2, y2 = points[-1]

    if abs(x1 - x2) < tol and abs(y1 - y2) < tol:
        return points[:-1]

    return points






# ═══════════════════════════════════════════════════════════════
# BLOC 2 — NOYAU GREEN ANALYTIQUE EXACT
# ═══════════════════════════════════════════════════════════════
#
# Théorème de Green :  ∬_Ω f dA = ∮_∂Ω Q(x,y) dy  avec ∂Q/∂x = f
#
# Arête paramétrée : x(s)=xa+s·Δx, y(s)=ya+s·Δy, s∈[0,1]
# ∮ Q dy = Σ_arêtes ∫₀¹ Q(x(s),y(s))·Δy ds
#
# ── PRIMITIVES Q (∂Q/∂x = f) ────────────────────────────────
#
# ══ ELS Zone C : f = C·ε = C·(ε₀+α·x+β·y) ══════════════════
#   Q_N  = C·(ε₀·x + α·x²/2 + β·x·y)
#   Q_My = C·(ε₀·x·y + α·x²·y/2 + β·x·y²)
#   Q_Mz = C·(ε₀·x²/2 + α·x³/3 + β·x²·y/2)
#   Q_Sc = x                                  [f=1 zone comprimée]
#
# Sur arête (s∈[0,1]) :
#   Notations : x=xa+s·Δx, y=ya+s·Δy
#   Moments : ∫s^k ds = 1/(k+1)
#   ∫x ds   = xa + Δx/2
#   ∫x² ds  = xa² + xa·Δx + Δx²/3
#   ∫x³ ds  = xa³ + 3xa²·Δx/2 + xa·Δx² + Δx³/4
#   ∫x·y ds = xa·ya + (xa·Δy+ya·Δx)/2 + Δx·Δy/3
#   ∫x²·y ds= xa²·ya + xa·(ya·Δx+xa·Δy)/1... (développé ci-dessous)
#   ∫x·y² ds= xa·ya² + (xa·2ya·Δy+ya²·Δx)/2 + xa·Δy²/3 + ya·Δx·Δy/2 + ...
#
#   N  = Δy · C · ∫₀¹(ε₀·x + α·x²/2 + β·x·y) ds
#      = Δy · C · [ε₀·(xa+Δx/2) + α/2·(xa²+xa·Δx+Δx²/3) + β·(xa·ya+(xa·Δy+ya·Δx)/2+Δx·Δy/3)]
#
#   My = Δy · C · ∫₀¹(ε₀·x·y + α·x²·y/2 + β·x·y²) ds
#
#   Mz = Δy · C · ∫₀¹(ε₀·x²/2 + α·x³/3 + β·x²·y/2) ds
#
# ══ ELU Zone R : f = fcd ════════════════════════════════════
#   Q_N  = fcd·x
#   Q_My = fcd·x·y
#   Q_Mz = fcd·x²/2
#   Q_Sc = x
#
#   N  = fcd · Δy · (xa + Δx/2)
#   My = fcd · Δy · (xa·ya + (xa·Δy+ya·Δx)/2 + Δx·Δy/3)
#   Mz = fcd · Δy/2 · (xa² + xa·Δx + Δx²/3)
#   Sc = Δy · (xa + Δx/2)
#
# ══ ELU Zone P : f = fcd·(2u-u²), u=ε/e2 ══════════════════
#   u = (ε₀+α·x+β·y)/e2 = u₀ + (α·x+β·y)/e2
#   u(s) = ua + s·Δu  où ua=(ε₀+α·xa+β·ya)/e2, Δu=α·Δx/e2+β·Δy/e2... 
#   NON : u dépend de x ET y simultanément via ε = ε₀+α·x+β·y
#
#   Sur l'arête : u(s) = (ε₀+α·(xa+s·Δx)+β·(ya+s·Δy))/e2
#               = ua + s·(α·Δx+β·Δy)/e2
#               = ua + s·Δu
#   C'est bien linéaire en s → u² est quadratique en s
#
#   f = fcd·(2u-u²) = fcd·(2(ua+s·Δu) - (ua+s·Δu)²)
#     = fcd·[(2ua-ua²) + s·(2Δu-2ua·Δu) - s²·Δu²]
#
#   Pour N = Δy·∫₀¹ Q_N(x(s),y(s)) ds, on a besoin de Q_N :
#   Q_N tel que ∂Q_N/∂x = fcd·(2u-u²)
#   Or u = (ε₀+α·x+β·y)/e2, donc ∂u/∂x = α/e2
#   ∂Q_N/∂x = fcd·(2u-u²)  →  Q_N = fcd·(u²/( α/e2) - u³/(3·α/e2))... 
#   c'est complexe. On utilise directement l'intégrale de ligne :
#
#   N  = ∫_arête σ·dy = Δy · ∫₀¹ fcd·(2ua+2Δu·s - ua²-2ua·Δu·s - Δu²·s²) ds
#      = Δy · fcd · [(2ua-ua²) + (2Δu-2ua·Δu)/2 - Δu²/3]
#      = Δy · fcd · [2ua - ua² + Δu - ua·Δu - Δu²/3]
#      = Δy · fcd · I0
#
#   My = ∫_arête σ·y·dy = Δy · ∫₀¹ fcd·(2u-u²)·(ya+Δy·s) ds
#      = Δy · fcd · [ya·I0 + Δy·I1]
#   où I1 = ∫₀¹(2u-u²)·s ds
#          = (2ua-ua²)/2 + (2Δu-2ua·Δu)/3 - Δu²/4
#
#   Mz = ∫_arête σ·x·dy = Δy · ∫₀¹ fcd·(2u-u²)·(xa+Δx·s) ds
#      = Δy · fcd · [xa·I0 + Δx·I1]
#
# IMPORTANT : Green donne ∬σ dA = ∮ Q dy
# Pour N : ∬σ dA = ∮ Q_N dy  NON IDENTIQUE à ∮ σ·dy  !!
# La formule correcte est :
#   ∬ f(x,y) dA = ∮_contour Q(x,y)·dy  avec ∂Q/∂x = f
#
# Pour ELU zone P, f=σ(x,y) et Q_N tel que ∂Q_N/∂x = σ(x,y) :
#   σ = fcd·(2ε/e2 - ε²/e2²)
#   ε = ε₀ + α·x + β·y
#   ∂σ/∂x = fcd·(2α/e2 - 2α·ε/e2²)
#   ∂Q_N/∂x = σ  →  Q_N = fcd·(2ε·x/(e2... NON
#   Q_N = ∫σ dx = fcd·∫(2(ε₀+αx+βy)/e2 - (ε₀+αx+βy)²/e2²) dx
#
#   Posons A = ε₀+βy (constant en x), B = α
#   ε = A + B·x
#   Q_N = fcd·∫(2(A+Bx)/e2 - (A+Bx)²/e2²) dx
#       = fcd·[2(Ax+Bx²/2)/e2 - (A²x+ABx²+B²x³/3)/e2²]
#       = fcd·[2Ax/e2 + Bx²/e2 - A²x/e2² - ABx²/e2² - B²x³/(3e2²)]
#
#   Sur l'arête, x=xa+s·Δx, y=ya+s·Δy :
#   A(s) = ε₀ + β·(ya+s·Δy)  → dépend de s via y !
#   Donc Q_N(x(s),y(s)) n'est pas simplement évaluable.
#
#   SOLUTION : On intègre directement ∫₀¹ Q_N(x(s),y(s))·Δy ds
#   en substituant x=xa+s·Δx et y=ya+s·Δy :
#   A(s) = ε₀+β·ya + β·Δy·s = A0 + dA·s  avec A0=ε₀+β·ya, dA=β·Δy
#   B = α (constant)
#   x(s) = xa + Δx·s
#
#   Q_N(s) = fcd·[2A(s)x(s)/e2 + B·x²(s)/e2
#                 - A²(s)·x(s)/e2² - A(s)B·x²(s)/e2² - B²·x³(s)/(3e2²)]
#
#   Chaque terme est un polynôme en s de degré ≤ 3 → intégrable exactement.
# ═══════════════════════════════════════════════════════════════

def _aire_signee(pts):
    """Aire signée — Shoelace. + = CCW, - = CW."""
    a = 0.0
    n = len(pts)
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i+1) % n]
        a += x1*y2 - x2*y1
    return 0.5*a


def _decouper_arete(xa, ya, xb, yb, eps0, alpha, beta, seuils):
    """
    Découpe l'arête aux t ∈ (0,1) où ε = seuil.
    Sécurisé pour les sections en I : préserve TOUJOURS l'orientation
    géométrique originale de l'arête pour ne pas fausser S_com.
    """
    ea = eps0 + alpha * xa + beta * ya
    eb = eps0 + alpha * xb + beta * yb
    de = eb - ea

    # On commence avec les bornes naturelles du segment dans le sens d'origine
    cuts = [0.0, 1.0]
    
    for s in seuils:
        if abs(de) > 1e-15:
            t = (s - ea) / de
            if 1e-11 < t < 1.0 - 1e-11:
                cuts.append(t)
                
    # TRÈS IMPORTANT : On trie TOUJOURS dans le sens croissant du paramètre t (0 -> 1)
    # afin de conserver le vecteur directeur géométrique original (xa,ya) -> (xb,yb)
    cuts = sorted(list(set(cuts)))

    segments = []
    for k in range(len(cuts) - 1):
        t1, t2 = cuts[k], cuts[k+1]
        tm = 0.5 * (t1 + t2)
        
        # Interpolation des coordonnées respectant le sens de parcours original
        x1 = xa + t1 * (xb - xa)
        y1 = ya + t1 * (yb - ya)
        x2 = xa + t2 * (xb - xa)
        y2 = ya + t2 * (yb - ya)
        
        # Déformation au milieu du sous-segment
        em = ea + tm * de
        
        segments.append((x1, y1, x2, y2, em))
        
    return segments
# ── Moments polynomiaux sur [0,1] ───────────────────────────────────
# On note p(s) = p0 + s·dp  pour tout paramètre linéaire en s.
# ∫₀¹ p^k · s^j ds se calcule par binôme.
#
# Fonctions utilitaires :
def _I(a, da, b, db, p, q):
    """
    ∫₀¹ (a+da·s)^p · (b+db·s)^q ds  pour p,q ∈ {0,1,2,3} entiers.
    Développé par binôme et intégré terme à terme.
    Retourne valeur scalaire.
    """
    # Développement binomial de (a+da·s)^p et (b+db·s)^q
    from math import comb
    total = 0.0
    for i in range(p+1):
        ci = comb(p, i) * (a**(p-i)) * (da**i)
        for j in range(q+1):
            cj = comb(q, j) * (b**(q-j)) * (db**j)
            total += ci * cj / (i+j+1)
    return total


def _contrib_Scom(x1, y1, x2, y2):
    """
    Calcul de l'aire comprimée par intégration directe (Gauss / Triangle).
    Marche parfaitement même sur contour ouvert (arêtes découpées).
    """
    return 0.5 * (x1 * y2 - x2 * y1)

def _contrib_ELS(x1, y1, x2, y2, eps0, alpha, beta, C):
    """
    Contribution ELS par intégration directe de la contrainte.
    Évite les erreurs de Green sur les contours ouverts après découpe.
    Convention : My = σ * y  et  Mz = σ * x
    """
    dx = x2 - x1
    dy = y2 - y1
    
    # Calcul des intégrales de base du triangle formé avec l'origine (0,0)
    # dA = 0.5 * (x1*y2 - x2*y1) -> On utilise l'approche de projection directe
    # Pour s'aligner sur ton noyau sans réécrire l'intégrateur, on conserve le dy de Green
    # mais on s'assure que l'arête respecte la fermeture.
    
    if abs(dy) < 1e-15:
        return np.zeros(3)

    Ix    = _I(x1, dx, 1, 0, 1, 0)   
    Ix2   = _I(x1, dx, 1, 0, 2, 0)   
    Ix3   = _I(x1, dx, 1, 0, 3, 0)   
    Ixy   = _I(x1, dx, y1, dy, 1, 1) 
    Ix2y  = _I(x1, dx, y1, dy, 2, 1) 
    Ixy2  = _I(x1, dx, y1, dy, 1, 2) 

    N_  = C * dy * (eps0 * Ix + (alpha / 2.0) * Ix2 + beta * Ixy)
    My_ = C * dy * (eps0 * Ixy + (alpha / 2.0) * Ix2y + beta * Ixy2)
    Mz_ = C * dy * ((eps0 / 2.0) * Ix2 + (alpha / 3.0) * Ix3 + (beta / 2.0) * Ix2y)

    return np.array([N_, My_, Mz_])





def _contrib_ELU_R(x1, y1, x2, y2, fcd):
    """
    Contribution EXACTE ELU zone R (σ=fcd) à [N, My, Mz].
    Q_N  = fcd·x        →  N  = fcd·Δy·∫x ds  = fcd·Δy·(x1+Δx/2)
    Q_My = fcd·x·y      →  My = fcd·Δy·∫x·y ds
    Q_Mz = fcd·x²/2     →  Mz = fcd·Δy/2·∫x² ds
    """
    dx = x2-x1;  dy = y2-y1
    if abs(dy) < 1e-15:
        return np.zeros(3)

    Ix  = _I(x1, dx, 1, 0, 1, 0)   # ∫x ds
    Ix2 = _I(x1, dx, 1, 0, 2, 0)   # ∫x² ds
    Ixy = _I(x1, dx, y1, dy, 1, 1) # ∫x·y ds

    N_  = fcd * dy * Ix
    My_ = fcd * dy * Ixy
    Mz_ = fcd * dy/2.0 * Ix2

    return np.array([N_, My_, Mz_])


def _contrib_ELU_P(x1, y1, x2, y2, eps0, alpha, beta, fcd, e2, n=2.0):
    """
    Contribution EXACTE ELU zone P (σ=fcd·(1 - (1-u)ⁿ), u=ε/e2) à [N, My, Mz].
    Généralisée pour une puissance n quelconque (ex: n=1.4 pour C70).
    
    Principe mathématique :
    Au lieu de développer le polynôme (impossible si n n'est pas entier),
    on intègre analytiquement la primitive exacte de la fonction puissance W^n.
    Pour éviter la division par zéro si alpha=0 ou beta=0, on projette
    le théorème de Green sur l'axe X ou Y selon le gradient dominant.
    """
    dx = x2 - x1
    dy = y2 - y1
    if abs(dx) < 1e-15 and abs(dy) < 1e-15:
        return np.zeros(3)

    # 1. Calcul des déformations aux nœuds de l'arête
    eps1 = eps0 + alpha * x1 + beta * y1
    eps2 = eps0 + alpha * x2 + beta * y2
    
    # 2. Variable interne W = 1 - u = 1 - eps / e2
    # Borne max à 0.0 pour éviter les NaN (numérique) avec des puissances non entières
    W1 = 1.0 - eps1 / e2
    W2 = 1.0 - eps2 / e2
    # Ne pas tronquer ici — poly_p est déjà clipé à eps <= e2
    # donc W1, W2 >= 0 par construction. La troncature masque des erreurs de clipping.
    # Si malgré tout W < 0 par tolérance numérique :
    W1 = max(W1, 0.0)
    W2 = max(W2, 0.0)
    # ← garder, MAIS s'assurer que la branche exacte est utilisée
    # quand dW est significatif, même si W1 = 0
    
    # Fonction d'intégration exacte de W(s)^p * s^k sur le segment paramétrique s ∈ [0, 1]
    def _I_W(p):
        dW = W2 - W1
        # Branche exacte prioritaire si dW est suffisant
        # ET si aucun des deux W n'est nul (risque NaN avec puissance négative en Taylor)
        use_exact = abs(dW) > 1e-3 or W1 == 0.0
        if not use_exact:
            # Taylor sécurisé — W1 > 0 garanti ici
            T0 = W1**p
            T1 = p * W1**(p-1) * dW
            T2 = p*(p-1)/2.0 * W1**(p-2) * dW**2
            T3 = p*(p-1)*(p-2)/6.0 * W1**(p-3) * dW**3
            I0 = T0 + T1/2 + T2/3 + T3/4
            I1 = T0/2 + T1/3 + T2/4 + T3/5
            return I0, I1
        else:
            # Formule exacte — valide même si W1 = 0
            I0 = (W2**(p+1) - W1**(p+1)) / (dW * (p+1))
            I1 = (W2**(p+2) - W1**(p+2)) / (dW**2 * (p+2)) - (W1/dW) * I0
            return I0, I1


    # 3. Cas particulier : Déformation uniforme sur toute la section (Compression pure)
    if abs(alpha) < 1e-12 and abs(beta) < 1e-12:
        sigma = fcd * (1.0 - W1**n)
        N_  = sigma * dy * (x1 + dx/2.0)
        My_ = sigma * dy * (x1*y1 + (x1*dy + y1*dx)/2.0 + dx*dy/3.0)
        Mz_ = sigma * dy * (x1**2 + x1*dx + dx**2/3.0) / 2.0
        return np.array([N_, My_, Mz_])

    # 4. Calcul des intégrales paramétriques I(p, k) nécessaires
    I0_n1, I1_n1 = _I_W(n + 1.0)
    I0_n2, _     = _I_W(n + 2.0)

    # 5. Application Adaptative du Théorème de Green
    # On choisit d'intégrer selon X ou Y pour ne JAMAIS diviser par un coefficient nul
    if abs(alpha) >= abs(beta):
        # --- Gradient dominant en X --- 
        # Vecteur de Green : Q_N = Int(sigma dx) => Contour sur dy
        C1 = e2 / (alpha * (n + 1.0))
        C2 = e2**2 / (alpha**2 * (n + 1.0) * (n + 2.0))
        
        # Parties purement géométriques/polynomiales
        P_N  = x1 + dx/2.0
        P_My = x1*y1 + (x1*dy + y1*dx)/2.0 + dx*dy/3.0
        P_Mz = (x1**2 + x1*dx + dx**2/3.0)/2.0
        
        # Assemblage avec les primitives exactes de la courbure
        N_  = fcd * dy * (P_N + C1 * I0_n1)
        My_ = fcd * dy * (P_My + C1 * (y1 * I0_n1 + dy * I1_n1))
        Mz_ = fcd * dy * (P_Mz + C1 * (x1 * I0_n1 + dx * I1_n1) + C2 * I0_n2)
        
    else:
        D1 = e2 / (beta * (n + 1.0))
        D2 = e2**2 / (beta**2 * (n + 1.0) * (n + 2.0))
    
        # Géométrie : en intégrant sur -dx, les rôles x↔y s'échangent
        P_N  = y1 + dy/2.0
        P_My = (y1**2 + y1*dy + dy**2/3.0) / 2.0   # était P_Mz en branche X
        P_Mz = x1*y1 + (x1*dy + y1*dx)/2.0 + dx*dy/3.0  # était P_My en branche X
    
        N_  = -fcd * dx * (P_N  + D1 * I0_n1)
        # ✅ D2 dans Mz (pas My) — symétrique de C2 dans Mz en branche X
        My_ = -fcd * dx * (P_My + D1 * (y1 * I0_n1 + dy * I1_n1) + D2 * I0_n2)
        Mz_ = -fcd * dx * (P_Mz + D1 * (x1 * I0_n1 + dx * I1_n1))
            
    return np.array([N_, My_, Mz_])




def domaines_comprimes(contour, eps0, alpha, beta):

    def eps(x,y):
        return eps0 + alpha*x + beta*y

    pts = np.asarray(contour)
    n = len(pts)

    domaines = []
    poly = []

    for i in range(n):

        xa, ya = pts[i]
        xb, yb = pts[(i+1) % n]

        ea = eps(xa, ya)
        eb = eps(xb, yb)

        # entrée en compression
        if ea <= 0 and eb > 0:
            t = -ea / (eb - ea)
            xi = xa + t*(xb-xa)
            yi = ya + t*(yb-ya)
            poly = [(xi, yi), (xb, yb)]

        # sortie compression
        elif ea > 0 and eb <= 0:
            t = -ea / (eb - ea)
            xi = xa + t*(xb-xa)
            yi = ya + t*(yb-ya)
            poly.append((xi, yi))
            domaines.append(poly)
            poly = []

        # entièrement en compression
        elif ea > 0 and eb > 0:
            if not poly:
                poly = [(xa, ya)]
            poly.append((xb, yb))

        # cas traction → rien

    # fermer si nécessaire
    if poly:
        domaines.append(poly)

    return domaines

from shapely.geometry import Polygon, LineString, MultiPolygon, GeometryCollection
from shapely.ops import split

def _orienter_polygone_ccw(pts):
    """
    Force un polygone orienté anti-horaire (CCW)
    """

    pts = np.asarray(pts)

    A = 0.0
    n = len(pts)

    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i+1) % n]
        A += x1*y2 - x2*y1

    if A < 0:
        return pts[::-1]   # inversion
    else:
        return pts

def _zones_comprimees(vertices, eps0, alpha, beta):

    poly = Polygon(vertices)

    if not poly.is_valid or poly.area == 0:
        return []

    # droite ε = 0
    L = 1e4

    if abs(beta) > 1e-12:
        x_vals = [-L, L]
        y_vals = [(-eps0 - alpha*x) / beta for x in x_vals]
    else:
        if abs(alpha) < 1e-12:
            return [poly] if eps0 > 0 else []
        x = -eps0 / alpha
        x_vals = [x, x]
        y_vals = [-L, L]

    line = LineString(list(zip(x_vals, y_vals)))

    try:
        result = split(poly, line)
    except:
        return []

    # ---------------------------------
    # ✅ gestion robuste des types
    # ---------------------------------
    zones = []

    def traiter_geom(g):
        if g.is_empty:
            return

        if g.geom_type == 'Polygon':
            cx, cy = g.centroid.x, g.centroid.y
            if eps0 + alpha*cx + beta*cy > 0:
                zones.append(g)

        elif g.geom_type in ['MultiPolygon', 'GeometryCollection']:
            for sub in g.geoms:
                traiter_geom(sub)

    traiter_geom(result)

    return zones

# =============================================================================
# DÉCOUPE GÉOMÉTRIQUE PAR PLANS DE DÉFORMATION
# =============================================================================

def _couper_poly_par_plan(poly, eps0, alpha, beta, valeur_seuil):
    """Coupe un polygone Shapely selon une ligne d'iso-déformation."""
    if not poly.is_valid or poly.area == 0:
        return []
    
    L = 1e4
    if abs(beta) > 1e-12:
        x_vals = [-L, L]
        y_vals = [(-eps0 - alpha*x + valeur_seuil) / beta for x in x_vals]
    else:
        if abs(alpha) < 1e-12:
            return [poly] if (eps0 > valeur_seuil) else []
        x = (valeur_seuil - eps0) / alpha
        x_vals, y_vals = [x, x], [-L, L]

    line = LineString(list(zip(x_vals, y_vals)))
    try:
        result = split(poly, line)
    except:
        return [poly]

    zones = []
    if result.geom_type == 'Polygon':
        zones.append(result)
    else:
        for g in result.geoms:
            if g.geom_type == 'Polygon' and not g.is_empty:
                zones.append(g)
    return zones

# =============================================================================
# MOTEUR D'INTÉGRATION PRINCIPAL
# =============================================================================
def clip_polygon_eps(vertices, eps0, alpha, beta, seuil, keep_above=True, tol=1e-12):
    """
    Coupe un polygone par un demi-plan défini par :
        eps(x,y) >= seuil   (si keep_above=True)
        eps(x,y) <= seuil   (si keep_above=False)

    Paramètres
    ----------
    vertices : array-like (N,2)
        Sommets du polygone (ordre quelconque mais fermé implicitement)
    eps0, alpha, beta : float
        Champ affine eps(x,y) = eps0 + alpha*x + beta*y
    seuil : float
        Valeur seuil
    keep_above : bool
        True  -> garde eps >= seuil
        False -> garde eps <= seuil
    tol : float
        Tolérance numérique

    Retour
    ------
    np.ndarray (M,2)
        Polygone clipé
    """

    def eps(x, y):
        return eps0 + alpha*x + beta*y

    def is_inside(e):
        if keep_above:
            return e >= seuil - tol
        else:
            return e <= seuil + tol

    output = []
    vertices = np.asarray(vertices, dtype=float)
    n = len(vertices)

    for i in range(n):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % n]

        e1 = eps(x1, y1)
        e2 = eps(x2, y2)

        inside1 = is_inside(e1)
        inside2 = is_inside(e2)

        # ✅ cas 1 : les deux points dedans
        if inside1 and inside2:
            output.append((x2, y2))

        # ✅ cas 2 : sortie du domaine
        elif inside1 and not inside2:
            if abs(e2 - e1) > tol:
                t = (seuil - e1) / (e2 - e1)
                xi = x1 + t * (x2 - x1)
                yi = y1 + t * (y2 - y1)
                output.append((xi, yi))

        # ✅ cas 3 : entrée dans le domaine
        elif not inside1 and inside2:
            if abs(e2 - e1) > tol:
                t = (seuil - e1) / (e2 - e1)
                xi = x1 + t * (x2 - x1)
                yi = y1 + t * (y2 - y1)
                output.append((xi, yi))
            output.append((x2, y2))

        # ✅ cas 4 : dehors → rien

    if len(output) == 0:
        return np.zeros((0, 2))

    return np.array(output)


def _integrer_polygone(contour, eps0, alpha, beta, mode, C=None, fck=None, fcd=None, e2=None):

    pts = _orienter_polygone_ccw(np.asarray(contour, float))

    # ✅ zone comprimée
    poly_c = clip_polygon_eps(pts, eps0, alpha, beta, 0.0)
    
    if len(poly_c) < 3:
        return np.zeros(4)

    res = np.zeros(4)

    if mode == 'ELS':

        edges = np.column_stack([poly_c, np.roll(poly_c, -1, axis=0)])

        for xa, ya, xb, yb in edges:
            dx = xb - xa
            dy = yb - ya

            res[3] += dy * (xa + dx/2)
            res[:3] += _contrib_ELS(xa, ya, xb, yb, eps0, alpha, beta, C)

        return res



    else:
        n= eps_n(fck)
        # ✅ zone rectangle
        poly_r = clip_polygon_eps(poly_c, eps0, alpha, beta, e2, keep_above=True)

        # ✅ zone parabole = poly_c - poly_r
        res_r = np.zeros(4)
        if len(poly_r) >= 3:
            edges = np.column_stack([poly_r, np.roll(poly_r, -1, axis=0)])
            for xa, ya, xb, yb in edges:
                dx = xb - xa
                dy = yb - ya
                res_r[3] += dy * (xa + dx/2)
                res_r[:3] += _contrib_ELU_R(xa, ya, xb, yb, fcd)

        # ✅ zone parabole (reste)
        poly_p = clip_polygon_eps(poly_c, eps0, alpha, beta, e2, keep_above=False)
        res_p = np.zeros(4)
        if len(poly_p) < 3:
            return res_r
        edges = np.column_stack([poly_p, np.roll(poly_p, -1, axis=0)])

        for xa, ya, xb, yb in edges:

            xm = 0.5*(xa+xb)
            ym = 0.5*(ya+yb)
            if eps0 + alpha*xm + beta*ym <= e2:

                dx = xb - xa
                dy = yb - ya
                res_p[3] += dy * (xa + dx/2)
                res_p[:3] += _contrib_ELU_P(xa, ya, xb, yb,
                                            eps0, alpha, beta, fcd, e2, n=n)


        return res_p + res_r



# ═══════════════════════════════════════════════════════════════
# BLOC 3 — FONCTIONS PUBLIQUES
# ═══════════════════════════════════════════════════════════════

@func
def S_com(polygon1, evi, eps0, alpha, beta):

    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

    res_c = _integrer_polygone(
        contour_cg, eps0, alpha, beta,
        mode='ELS', C=1.0
    )

    res_v = np.zeros(4)

    for trou in evidements_cg:
        res_v += _integrer_polygone(
            trou, eps0, alpha, beta,
            mode='ELS', C=1.0
        )

    return res_c[3] - res_v[3]



@func
def calculer_N_My_Mz(
    polygon1, evi, p_acier, s_acier,
    a_com, n, eps0, alpha, beta
):
    """ELS — N, My, Mz par Green analytique EXACT. Sans nettoyer_polygone."""
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return 0.0, 0.0, 0.0

    C = 200000.0 / (float(n) * 1000.0)

    res_c = _integrer_polygone(contour_cg, eps0, alpha, beta, 'ELS', C=C)
    res_v = np.zeros(4)
    for trou in evidements_cg:
        if len(trou) >= 3:
            res_v += _integrer_polygone(trou, eps0, alpha, beta, 'ELS', C=C)

    Nc, Mc_y, Mc_z = (res_c - res_v)[:3]

    Ns = Msy = Msz = 0.0
    if p_acier is not None and len(p_acier) > 0:
        acier_data = acier_G(p_acier, s_acier)
        if len(acier_data) > 0:
            acier = np.array(acier_data, dtype=float)
            x_s = acier[:,0] - Cx;  y_s = acier[:,1] - Cy
            eps_s = eps0 + alpha*x_s + beta*y_s
            sig_s = sigma_s_lin1(eps_s, a_com)
            forces = sig_s * acier[:,2]
            Ns  = np.sum(forces)
            Msy = np.sum(forces * y_s)
            Msz = np.sum(forces * x_s)

    return Ns+Nc, Msy+Mc_y, Msz+Mc_z


@func
def solve_GG_ELS(polygon1, evi, p_acier, s_acier,
                 a_com, n, Nobj, Myobj, Mzobj):
    x0 = np.array([0.0, 0.001, 0.0])
    def residuals(eps):
        N, My, Mz = calculer_N_My_Mz(
            polygon1, evi, p_acier, s_acier, a_com, n, *eps)
        return np.array([N-Nobj, My-Myobj, Mz-Mzobj])
    def jacobian(eps, h=1e-6):
        J = np.empty((3,3)); r0 = residuals(eps)
        for i in range(3):
            dh = np.zeros(3); dh[i] = h
            J[:,i] = (residuals(eps+dh) - r0) / h
        return J
    result = root(residuals, x0, jac=jacobian, method='hybr',
                  tol=1e-6, options={'maxfev': 1000})
    return result.x


@func
def resultats_GG_ELS(
    polygon1, evi, p_acier, s_acier,
    a_com, n, eps0, alpha, beta
):
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return {}

    poly_array = np.asarray(contour_cg, dtype=float)
    eps_c = eps0 + alpha*poly_array[:,0] + beta*poly_array[:,1]
    espcmax = float(np.max(eps_c)); espcmin = float(np.min(eps_c))
    sig_c_max = float(sigma_c_n1(espcmax, n))
    sig_c_min = float(sigma_c_n1(espcmin, n))

    espsmax = espsmin = sig_s_max = sig_s_min = ""
    if p_acier is not None and len(p_acier) > 0:
        acier = acier_G(p_acier, s_acier)
        if len(acier) > 0:
            acier_array = np.asarray(acier, dtype=float)
            x_s = acier_array[:,0] - Cx; y_s = acier_array[:,1] - Cy
            eps_s = eps0 + alpha*x_s + beta*y_s
            espsmax = float(np.max(eps_s)); espsmin = float(np.min(eps_s))
            sig_s_max = float(sigma_s_lin1(espsmax, a_com))
            sig_s_min = float(sigma_s_lin1(espsmin, a_com))

    A_com = S_com(polygon1, evi, eps0, alpha, beta)
    N, My, Mz = calculer_N_My_Mz(
        polygon1, evi, p_acier, s_acier, a_com, n, eps0, alpha, beta)

    return {
        'ACOM': A_com,
        'EPS_C_MAX': espcmax,  'EPS_C_MIN': espcmin,
        'SIG_C_MAX': sig_c_max, 'SIG_C_MIN': sig_c_min,
        'EPS_S_MAX': espsmax,  'EPS_S_MIN': espsmin,
        'SIG_S_MAX': sig_s_max, 'SIG_S_MIN': sig_s_min,
        'N': N, 'MY': My, 'MZ': Mz, 'CX': Cx, 'CY': Cy
    }


@func
def e_resultats_GG_ELS(polygon1, evi, p_acier, s_acier,
                        a_com, n, eps0, alpha, beta, resultats):
    tout = resultats_GG_ELS(
        polygon1, evi, p_acier, s_acier, a_com, n, eps0, alpha, beta)
    return [tout[r] for r in resultats.split(',')]


@func
def calculer_N_My_Mz_ELU_pararect(
    polygon1, evi, p_acier, s_acier,
    a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
    eps0, alpha, beta
):
    """ELU — N, My, Mz par Green analytique EXACT. Sans nettoyer_polygone."""
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return 0.0, 0.0, 0.0
    n = float(eps_n(fck))
    e2 = float(eps_c2(fck))

    res_c = _integrer_polygone(contour_cg, eps0, alpha, beta,
                               'ELU',fck=fck, fcd=fcd, e2=e2)
    res_v = np.zeros(4)
    for trou in evidements_cg:
        if len(trou) >= 3:
            res_v += _integrer_polygone(trou, eps0, alpha, beta,
                                        'ELU',fck=fck, fcd=fcd, e2=e2)

    Nc, Mc_y, Mc_z = (res_c - res_v)[:3]

    Ns = Msy = Msz = 0.0
    if p_acier is not None and len(p_acier) > 0:
        acier_data = acier_G(p_acier, s_acier)
        if len(acier_data) > 0:
            acier = np.array(acier_data, dtype=float)
            x_s = acier[:,0] - Cx;  y_s = acier[:,1] - Cy
            eps_s = eps0 + alpha*x_s + beta*y_s
            sig_s = sigma_s_palier1(fyd, k, eps_uk, eps_ud, eps_s, a_com)
            forces = sig_s * acier[:,2]
            Ns  = np.sum(forces)
            Msy = np.sum(forces * y_s)
            Msz = np.sum(forces * x_s)

    return Ns+Nc, Msy+Mc_y, Msz+Mc_z


@func
def solve_GG_ELU_pararect(polygon1, evi, p_acier, s_acier,
                           a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
                           Nobj, Myobj, Mzobj):
    x0 = np.array([0.0, 0.1, 0.0])
    def residuals(eps):
        N, My, Mz = calculer_N_My_Mz_ELU_pararect(
            polygon1, evi, p_acier, s_acier,
            a_com, fck, fcd, fyd, k, eps_uk, eps_ud, *eps)
        return np.array([N-Nobj, My-Myobj, Mz-Mzobj])
    def jacobian(eps, h=1e-6):
        J = np.empty((3,3)); r0 = residuals(eps)
        for i in range(3):
            dh = np.zeros(3); dh[i] = h
            J[:,i] = (residuals(eps+dh) - r0) / h
        return J
    #result = root(residuals, x0, jac=jacobian, method='hybr',
    #              tol=1e-6, options={'maxfev': 1000})
    #result = root(residuals, x0, method='df-sane',
    #          tol=1e-6, options={'maxfev': 2000, 'fatol': 1e-6})
    #from scipy.optimize import fsolve
    sol, info, ier, msg = fsolve(residuals, x0, full_output=True, 
                              xtol=1e-6)
    #result = root(residuals, x0, method='krylov',
    #          tol=1e-5, options={'maxiter': 200})
    return sol #result.x




@func
def resultats_GG_ELU_pararect(
    polygon1, evi, p_acier, s_acier,
    a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
    eps0, alpha, beta
):
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return {}

    poly_array = np.asarray(contour_cg, dtype=float)
    eps_c = eps0 + alpha*poly_array[:,0] + beta*poly_array[:,1]
    espcmax = float(np.max(eps_c)); espcmin = float(np.min(eps_c))
    sig_c_max = float(sigma_c_pararect1(fck, fcd, espcmax))
    sig_c_min = float(sigma_c_pararect1(fck, fcd, espcmin))

    epssmax = epssmin = sig_s_max = sig_s_min = ""
    pour_a_pl = 0.0
    if p_acier is not None and len(p_acier) > 0:
        acier = acier_G(p_acier, s_acier)
        if len(acier) > 0:
            acier_array = np.asarray(acier, dtype=float)
            x_s = acier_array[:,0] - Cx; y_s = acier_array[:,1] - Cy
            eps_s = eps0 + alpha*x_s + beta*y_s
            epssmax = float(np.max(eps_s)); epssmin = float(np.min(eps_s))
            plast = np.sum(np.abs(eps_s) >= (float(fyd)/200.0))
            pour_a_pl = float(plast/len(eps_s)*100.0)
            sig_s_max = float(sigma_s_palier1(fyd, k, eps_uk, eps_ud, epssmax, a_com))
            sig_s_min = float(sigma_s_palier1(fyd, k, eps_uk, eps_ud, epssmin, a_com))

    A_com = float(S_com(polygon1, evi, eps0, alpha, beta))
    N, My, Mz = calculer_N_My_Mz_ELU_pararect(
        polygon1, evi, p_acier, s_acier,
        a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
        eps0, alpha, beta)

    return {
        'EPS0': eps0,   'ALPHA': alpha,   'BETA': beta,
        'ACOM': A_com,
        'EPS_C_MAX': espcmax,   'EPS_C_MIN': espcmin,
        'SIG_C_MAX': sig_c_max, 'SIG_C_MIN': sig_c_min,
        'EPS_S_MAX': epssmax,   'EPS_S_MIN': epssmin,
        'SIG_S_MAX': sig_s_max, 'SIG_S_MIN': sig_s_min,
        'N': float(N),  'MY': float(My),  'MZ': float(Mz),
        'PA': pour_a_pl, 'CX': Cx,        'CY': Cy
    }


@func
def e_resultats_GG_ELU_pararect(polygon1, evi, p_acier, s_acier,
                                  a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
                                  eps0, alpha, beta, resultats):
    tout = resultats_GG_ELU_pararect(
        polygon1, evi, p_acier, s_acier,
        a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
        eps0, alpha, beta)
    return [tout[r] for r in resultats.split(',')]


