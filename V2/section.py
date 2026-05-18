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




# section.py — version complète avec intégration analytique exacte
# par théorème de Green pour loi parabole-rectangle n=2 (C20-C50)

import numpy as np
from scipy.optimize import root, fsolve
from matplotlib.path import Path
from scipy.spatial import Delaunay

# ══════════════════════════════════════════════════════════════════════
# FONCTIONS DE BASE
# ══════════════════════════════════════════════════════════════════════

def sgn(x):
    return np.where(x < 0, -1, 1)

def eps_c2(fck_array):
    """Renvoie la déformation eps_c2 selon NF EN 1992-1-1 Tab. 3.1, en mm/m"""
    fck = np.asarray(fck_array)
    return np.where(fck <= 50, 2.0, 2.0 + 0.085 * (fck - 50)**0.53)

def eps_n(fck_array):
    """Renvoie le coefficient n selon NF EN 1992-1-1 Tab. 3.1"""
    fck = np.asarray(fck_array)
    return np.where(fck <= 50, 2.0, 1.4 + 23.4 * ((90 - fck) / 100)**4)

def sigma_c_n(eps_c_array, n_array):
    """Contrainte béton fissuré ELS — loi linéaire."""
    eps   = np.asarray(eps_c_array)
    n     = np.asarray(n_array)
    sigma = (200000 / n) * (eps / 1000)
    return np.maximum(0, sigma)

def eps_c_n(sig_c_array, n_array):
    """Déformation béton fissuré ELS — inverse loi linéaire."""
    sig     = np.asarray(sig_c_array)
    n       = np.asarray(n_array)
    epsilon = (1000 * n / 200000) * sig
    return np.where(sig > 0, epsilon, 0.0)

def sigma_c_pararect(fck, fcd, eps_c_array):
    """Loi parabole-rectangle EC2 § 3.1.7 — vectorisée."""
    eps_c      = np.asarray(eps_c_array)
    e2         = eps_c2(fck)
    n          = eps_n(fck)
    conditions = [eps_c <= 0,
                  (eps_c > 0) & (eps_c <= e2),
                  eps_c > e2]
    choix      = [0.0,
                  fcd * (1 - (1 - eps_c / e2)**n),
                  fcd]
    return np.select(conditions, choix)

def epsilon_c_pararect(fck, fcd, sig_c_array):
    """Inverse loi parabole-rectangle."""
    sig      = np.asarray(sig_c_array)
    e2       = eps_c2(fck)
    n        = eps_n(fck)
    sig_safe = np.clip(sig, 0, fcd * 0.9999999)
    epsilon  = e2 * (1 - (1 - sig_safe / fcd)**(1/n))
    return np.select([sig <= 0, sig >= fcd], [0.0, e2], default=epsilon)

def sigma_s_palier(fyd, k, eps_uk, eps_ud, eps_s_array):
    """Contrainte acier — loi bilinéaire avec palier ELU."""
    eps_s  = np.array(eps_s_array)
    eps_sd = fyd / 200000 * 1000
    sigma1 = 200000 * eps_s / 1000
    sigma2 = sgn(eps_s) * (fyd + (k * fyd - fyd) / (eps_uk - eps_sd) * (np.abs(eps_s) - eps_sd))
    return np.where(np.abs(eps_s) > eps_sd, sigma2, sigma1)

def sigma_s_palier1(fyd, k, eps_uk, eps_ud, eps_s_array, a_com):
    """Contrainte acier avec pondération compression ELU."""
    eps_s  = np.array(eps_s_array)
    sigma1 = sigma_s_palier(fyd, k, eps_uk, eps_ud, eps_s) * a_com
    sigma2 = sigma_s_palier(fyd, k, eps_uk, eps_ud, eps_s)
    return np.where(eps_s > 0, sigma1, sigma2)

def eps_s_palier(fyd, k, eps_uk, eps_ud, sig_s_array):
    """Inverse loi bilinéaire acier."""
    sig      = np.asarray(sig_s_array)
    Es       = 200000
    eps_sd   = (fyd / Es) * 1000
    pente_inv = (eps_uk - eps_sd) / (k * fyd - fyd) if k > 1 else 0
    abs_sig   = np.abs(sig)
    eps_elaste = (sig / Es) * 1000
    eps_plast  = np.sign(sig) * (eps_sd + (abs_sig - fyd) * pente_inv)
    conditions = [(sig == 0), (abs_sig <= fyd),
                  (abs_sig > fyd) & (abs_sig <= k * fyd)]
    return np.select(conditions, [0.0, eps_elaste, eps_plast], default=0.0)

def sigma_s_lin(eps_s_array):
    """Contrainte acier — loi élastique linéaire infinie ELS."""
    return 200000 * np.asarray(eps_s_array) / 1000

def eps_s_lin(sig_s_array):
    """Déformation acier — inverse loi linéaire ELS."""
    return 1000 * np.asarray(sig_s_array) / 200000

def sigma_s_lin1(eps_s, a_com):
    """Contrainte acier linéaire avec pondération compression."""
    eps      = np.asarray(eps_s)
    sig_base = 200000 * eps / 1000
    return np.where(eps > 0.0, sig_base * a_com, sig_base)

def sigma_c_n1(eps_c, n):
    """Contrainte béton fissuré sans np.where."""
    eps     = np.asarray(eps_c)
    eps_abs = np.abs(eps)
    filtre  = (eps + eps_abs) / (2 * eps_abs + 1e-30)
    return (200000 / n * eps / 1000) * filtre

def sigma_c_pararect1(fck, fcd, eps_c_array):
    """Variante loi parabole-rectangle sans np.select."""
    eps_c      = np.array(eps_c_array)
    eps_c2_val = eps_c2(fck)
    eps_n_val  = eps_n(fck)
    term1 = np.abs(1 - eps_c / eps_c2_val)
    term2 = (1 - eps_c / eps_c2_val)
    sigma = (fcd * (1 - (term1 + term2) / (2 * term1 + 1e-30) * term1**eps_n_val)
             * (eps_c + np.abs(eps_c)) / (2 * np.abs(eps_c) + 1e-30))
    return sigma


# ══════════════════════════════════════════════════════════════════════
# INTÉGRATION ANALYTIQUE EXACTE PAR THÉORÈME DE GREEN — n=2
# ══════════════════════════════════════════════════════════════════════
#
# Champ de déformation : ε(x,y) = ε₀ + α·x + β·y
#
# Loi béton n=2 :
#   σ(ε) = 0                              si ε ≤ 0
#   σ(ε) = fcd·(2ε/e2 - ε²/e2²)         si 0 < ε ≤ e2
#   σ(ε) = fcd                            si ε > e2
#
# Théorème de Green : ∬_Ω f dA = ∮_∂Ω H dy   avec ∂H/∂x = f
#
# Pour N = ∬σ dA  : on cherche H_N tel que ∂H_N/∂x = σ(ε₀+αx+βy)
# Pour My = ∬σy dA : H_My tel que ∂H_My/∂x = σ·y
# Pour Mz = ∬σx dA : H_Mz tel que ∂H_Mz/∂x = σ·x
#
# Sur une arête (x1,y1)→(x2,y2), on paramètre :
#   x(t) = x1 + t·dx,  y(t) = y1 + t·dy,  t∈[0,1]
#   ε(t) = ε_a + α·dx·t  où ε_a = ε₀ + α·x1 + β·y(t)
#
# L'intégrale de ligne ∫₀¹ H(x(t),y(t))·dy est calculée EXACTEMENT
# car H est un polynôme en x (degré 3 au max pour la zone parabolique).
#
# ── Primitives de H par rapport à t sur l'arête ──────────────────────
#
# Zone parabolique (0 < ε ≤ e2) :
#   σ = fcd·(2/e2·ε - 1/e2²·ε²)
#   ε = ε_a + α·dx·t
#
#   ∫σ·dy = fcd·dy·∫₀¹(2/e2·(ε_a+α·dx·t) - 1/e2²·(ε_a+α·dx·t)²) dt
#
# Intégrales élémentaires :
#   ∫₀¹ t^k dt = 1/(k+1)
#   ∫₀¹ (ε_a+c·t)   dt = ε_a + c/2
#   ∫₀¹ (ε_a+c·t)²  dt = ε_a² + ε_a·c + c²/3
#   ∫₀¹ (ε_a+c·t)³  dt = ε_a³ + 3ε_a²·c/2 + ε_a·c²+ c³/4  [pour Mz]
#
# Zone rectangulaire (ε > e2) :
#   σ = fcd  →  ∫σ·dy = fcd·dy·1 = fcd·dy
#   ∫σ·y·dy = fcd·dy·(y1+dy/2)    [y moyen sur l'arête]
#   ∫σ·x·dy = fcd·dy·(x1+dx/2)    [x moyen sur l'arête]
#
# MAIS : sur une arête, ε varie linéairement → les zones changent
# en cours d'arête. On décompose donc chaque arête en sous-intervalles
# [0,t_a], [t_a,t_b], [t_b,1] selon les racines ε=0 et ε=e2.
# Sur chaque sous-intervalle la zone est constante → intégrale EXACTE.
# ══════════════════════════════════════════════════════════════════════

def _racines_segment(eps_a, eps_b, seuils):
    """
    Trouve les paramètres t∈(0,1) où ε(t) = ε_a + (ε_b-ε_a)*t
    franchit les valeurs dans seuils.
    Retourne liste triée de t∈(0,1).
    """
    racines = []
    deps = eps_b - eps_a
    if abs(deps) < 1e-30:
        return racines
    for s in seuils:
        t = (s - eps_a) / deps
        if 1e-12 < t < 1 - 1e-12:
            racines.append(t)
    return sorted(racines)


def _integrale_exacte_zone(eps_a, deps, fcd, e2, dy, y1, dy_seg, x1, dx_seg, zone):
    """
    Intégrale EXACTE sur un sous-intervalle [0,1] complet
    (après reparamétrisation) pour N, My, Mz.

    ε(t) = eps_a + deps*t  sur [0,1]
    x(t) = x1 + dx_seg*t
    y(t) = y1 + dy_seg*t

    N  = ∫₀¹ σ(ε) · dy_seg dt
    My = ∫₀¹ σ(ε) · y(t) · dy_seg dt
    Mz = ∫₀¹ σ(ε) · x(t) · dy_seg dt

    zone : 0 = nulle, 1 = parabolique, 2 = rectangulaire
    """
    if zone == 0:
        return 0.0, 0.0, 0.0

    if zone == 2:
        # σ = fcd = constante
        # ∫₀¹ fcd · dy_seg dt = fcd · dy_seg
        N_  = fcd * dy_seg
        # ∫₀¹ fcd · (y1 + dy_seg·t) · dy_seg dt = fcd·dy_seg·(y1 + dy_seg/2)
        My_ = fcd * dy_seg * (y1 + dy_seg / 2.0)
        # ∫₀¹ fcd · (x1 + dx_seg·t) · dy_seg dt = fcd·dy_seg·(x1 + dx_seg/2)
        Mz_ = fcd * dy_seg * (x1 + dx_seg / 2.0)
        return N_, My_, Mz_

    # zone == 1 : σ = fcd·(2ε/e2 - ε²/e2²)  avec ε = eps_a + deps·t
    # ── Moments en ε : ──────────────────────────────────────────────
    # ∫₀¹ ε dt        = eps_a + deps/2
    # ∫₀¹ ε² dt       = eps_a² + eps_a·deps + deps²/3
    # ∫₀¹ ε·t dt      = eps_a/2 + deps/3
    # ∫₀¹ ε²·t dt     = eps_a²/2 + 2·eps_a·deps/3 + deps²/4  [pour Mz]

    I_e  = eps_a + deps / 2.0
    I_e2 = eps_a**2 + eps_a * deps + deps**2 / 3.0

    # N = fcd·dy_seg·∫₀¹ (2ε/e2 - ε²/e2²) dt
    N_ = fcd * dy_seg * (2.0 * I_e / e2 - I_e2 / e2**2)

    # My = fcd·dy_seg·∫₀¹ (2ε/e2 - ε²/e2²)·(y1 + dy_seg·t) dt
    # = fcd·dy_seg·[ (2/e2)·∫ε·(y1+dy_seg·t) - (1/e2²)·∫ε²·(y1+dy_seg·t) ]
    # ∫₀¹ ε·(y1+dy_seg·t) dt = y1·I_e + dy_seg·∫₀¹ε·t dt
    #                         = y1·I_e + dy_seg·(eps_a/2 + deps/3)
    # ∫₀¹ ε²·(y1+dy_seg·t) dt= y1·I_e2 + dy_seg·∫₀¹ε²·t dt
    #                         = y1·I_e2 + dy_seg·(eps_a²/2 + 2·eps_a·deps/3 + deps²/4)

    I_et  = eps_a / 2.0 + deps / 3.0
    I_e2t = eps_a**2 / 2.0 + 2.0 * eps_a * deps / 3.0 + deps**2 / 4.0

    int_e_y  = y1 * I_e  + dy_seg * I_et
    int_e2_y = y1 * I_e2 + dy_seg * I_e2t

    My_ = fcd * dy_seg * (2.0 * int_e_y / e2 - int_e2_y / e2**2)

    # Mz = fcd·dy_seg·∫₀¹ (2ε/e2 - ε²/e2²)·(x1 + dx_seg·t) dt
    # même structure avec x à la place de y
    int_e_x  = x1 * I_e  + dx_seg * I_et
    int_e2_x = x1 * I_e2 + dx_seg * I_e2t

    Mz_ = fcd * dy_seg * (2.0 * int_e_x / e2 - int_e2_x / e2**2)

    return N_, My_, Mz_


def _green_segment_exact(x1, y1, x2, y2, fcd, e2, eps0, alpha, beta):
    """
    Contribution EXACTE d'un segment au torseur [N, My, Mz]
    par théorème de Green — loi parabole-rectangle n=2.

    Décompose le segment aux points de transition ε=0 et ε=e2,
    puis intègre analytiquement sur chaque sous-segment.
    """
    dx = x2 - x1
    dy = y2 - y1
    if abs(dx) < 1e-15 and abs(dy) < 1e-15:
        return 0.0, 0.0, 0.0

    eps_debut = eps0 + alpha * x1 + beta * y1
    eps_fin   = eps0 + alpha * x2 + beta * y2

    # Points de transition t ∈ (0,1)
    t_transitions = [0.0] + _racines_segment(eps_debut, eps_fin, [0.0, e2]) + [1.0]

    N_tot = My_tot = Mz_tot = 0.0

    for k in range(len(t_transitions) - 1):
        ta = t_transitions[k]
        tb = t_transitions[k + 1]
        dt = tb - ta
        if dt < 1e-15:
            continue

        # Reparamétrisation sur [0,1] : τ → ta + dt·τ
        # ε(τ) = ε(ta) + (ε(tb)-ε(ta))·τ
        eps_a = eps_debut + (eps_fin - eps_debut) * ta
        eps_b = eps_debut + (eps_fin - eps_debut) * tb
        deps  = eps_b - eps_a

        # Coordonnées au début du sous-segment (en τ=0)
        x_a = x1 + ta * dx
        y_a = y1 + ta * dy
        # Incréments du sous-segment
        dx_s = dx * dt
        dy_s = dy * dt

        # Détermination de la zone (au milieu du sous-segment)
        eps_mid = (eps_a + eps_b) / 2.0
        if eps_mid <= 0:
            zone = 0
        elif eps_mid <= e2:
            zone = 1
        else:
            zone = 2

        N_, My_, Mz_ = _integrale_exacte_zone(
            eps_a, deps, fcd, e2,
            dy_s, y_a, dy_s, x_a, dx_s,
            zone
        )
        N_tot  += N_
        My_tot += My_
        Mz_tot += Mz_

    return N_tot, My_tot, Mz_tot


def _aire_signee_polygone(pts):
    """Aire signée (shoelace). Positif = CCW, négatif = CW."""
    pts = np.asarray(pts, dtype=float)
    n   = len(pts)
    a   = 0.0
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        a += x1 * y2 - x2 * y1
    return 0.5 * a


def _integrer_polygone_exact(polygon_pts, fcd, e2, eps0, alpha, beta):
    """
    Intégrale EXACTE de [N, My, Mz] sur un polygone par Green.
    Indépendant du sens de saisie — correction automatique du signe.
    """
    pts   = np.asarray(polygon_pts, dtype=float)
    n_pts = len(pts)
    if n_pts < 3:
        return 0.0, 0.0, 0.0

    N_tot = My_tot = Mz_tot = 0.0

    for i in range(n_pts):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n_pts]
        dN, dMy, dMz = _green_segment_exact(
            x1, y1, x2, y2, fcd, e2, eps0, alpha, beta
        )
        N_tot  += dN
        My_tot += dMy
        Mz_tot += dMz

    # Correction selon le sens de saisie
    aire = _aire_signee_polygone(pts)
    if aire < 0:
        N_tot  = -N_tot
        My_tot = -My_tot
        Mz_tot = -Mz_tot

    return N_tot, My_tot, Mz_tot


# ══════════════════════════════════════════════════════════════════════
# INTÉGRATION APPROCHÉE (Gauss-Delaunay) — conservée pour n ≠ 2
# ══════════════════════════════════════════════════════════════════════

def polygone_integrate(f, vertices, tol=1e-9, rtol=1e-9, max_depth=3):
    """
    Intégrateur adaptatif Gauss-Delaunay.
    Utilisé pour les lois non-polynomiales (n ≠ 2, ELS linéaire, etc.)
    """
    if not hasattr(polygone_integrate, '_gauss'):
        pts = np.array([
            [0.3333333333333333, 0.3333333333333333],
            [0.1012865073234099, 0.1012865073234099],
            [0.7974269853530873, 0.1012865073234099],
            [0.1012865073234099, 0.7974269853530873],
            [0.4701420641051151, 0.0597158717897698],
            [0.4701420641051151, 0.4701420641051151],
            [0.0597158717897698, 0.4701420641051151],
        ], dtype=np.float64)
        w = np.array([0.225,
                      0.1259391805448272, 0.1259391805448272, 0.1259391805448272,
                      0.1323941527885062, 0.1323941527885062, 0.1323941527885062],
                     dtype=np.float64)
        polygone_integrate._gauss = (pts.T, w)

    gpts_T, gw = polygone_integrate._gauss

    def gauss_integral(v0, v1, v2, area):
        J0x, J0y = v1[0]-v0[0], v1[1]-v0[1]
        J1x, J1y = v2[0]-v0[0], v2[1]-v0[1]
        X = v0[0] + J0x*gpts_T[0] + J1x*gpts_T[1]
        Y = v0[1] + J0y*gpts_T[0] + J1y*gpts_T[1]
        valeurs = np.asarray(f(X, Y))
        if valeurs.ndim > 1:
            return area * np.dot(gw, valeurs.T)
        return area * np.dot(gw, valeurs)

    def integrate_triangle_iter(v0, v1, v2, atol):
        stack = [(v0, v1, v2, atol, 0)]
        total_tri = None
        while stack:
            v0, v1, v2, atol, depth = stack.pop()
            area = 0.5 * abs((v1[0]-v0[0])*(v2[1]-v0[1])-(v2[0]-v0[0])*(v1[1]-v0[1]))
            if area < 1e-15:
                continue
            coarse = gauss_integral(v0, v1, v2, area)
            if total_tri is None:
                total_tri = np.zeros_like(coarse)
            if depth >= max_depth:
                total_tri += coarse
                continue
            m01, m12, m02 = 0.5*(v0+v1), 0.5*(v1+v2), 0.5*(v0+v2)
            a4 = area * 0.25
            g0 = gauss_integral(v0, m01, m02, a4)
            g1 = gauss_integral(v1, m01, m12, a4)
            g2 = gauss_integral(v2, m12, m02, a4)
            g3 = gauss_integral(m01, m12, m02, a4)
            fine = g0+g1+g2+g3
            err  = np.max(np.abs(fine - coarse))
            scale = max(np.max(np.abs(fine)), 1e-30)
            if err <= atol and err <= rtol * scale:
                total_tri += fine
            else:
                atol4 = atol * 0.25
                stack.extend([
                    (v0, m01, m02, atol4, depth+1),
                    (v1, m01, m12, atol4, depth+1),
                    (v2, m12, m02, atol4, depth+1),
                    (m01, m12, m02, atol4, depth+1)
                ])
        return total_tri if total_tri is not None else 0.0

    if vertices is None or len(vertices) < 3:
        return 0.0

    vertices     = np.asarray(vertices, dtype=np.float64)
    tri          = Delaunay(vertices)
    polygon_path = Path(vertices)
    total_final  = None

    for simplex in tri.simplices:
        v0, v1, v2 = vertices[simplex]
        if polygon_path.contains_point((v0+v1+v2)/3.0):
            res = integrate_triangle_iter(v0, v1, v2, tol)
            if total_final is None:
                total_final = np.zeros_like(res)
            total_final += res

    return total_final if total_final is not None else 0.0


# ══════════════════════════════════════════════════════════════════════
# GESTION DES ARMATURES
# ══════════════════════════════════════════════════════════════════════

def acier_G(p_acier, s_acier):
    """Construit la liste des armatures [x, y, section en cm²]."""
    if p_acier is None or len(p_acier) == 0:
        return []
    p_acier = np.atleast_2d(p_acier)
    s_acier = np.atleast_1d(s_acier)
    newacier = []
    nx = min(p_acier.shape[0], s_acier.shape[0])
    for i in range(nx):
        try:
            if p_acier[i, 0] in (None, "", " ") or p_acier[i, 1] in (None, "", " "):
                continue
            x = float(p_acier[i, 0])
            y = float(p_acier[i, 1])
            s = float(s_acier[i])
            newacier.append([x, y, s / 10000.0])
        except (ValueError, TypeError, IndexError):
            continue
    return newacier


# ══════════════════════════════════════════════════════════════════════
# GESTION DES LISTES ET POLYGONES
# ══════════════════════════════════════════════════════════════════════

def n_liste(data):
    """Compte le nombre de séparateurs 'fin'."""
    count = 0
    for row in data:
        if len(row) >= 2:
            v0 = str(row[0]).strip().lower()
            v1 = str(row[1]).strip().lower()
            if v0 == "fin" or v1 == "fin":
                count += 1
    return count

def liste_i(data, i):
    """Retourne la i-ème sous-liste séparée par 'fin'."""
    fins = []
    for idx, row in enumerate(data):
        if len(row) >= 2:
            v0 = str(row[0]).strip().lower()
            v1 = str(row[1]).strip().lower()
            if v0 == "fin" or v1 == "fin":
                fins.append(idx)
    if i < 1 or i > len(fins):
        return []
    start = 0 if i == 1 else fins[i-2] + 1
    end   = fins[i-1]
    return data[start:end]

def points_valides(data):
    """Filtre les points numériques valides."""
    pts = []
    for row in data:
        if len(row) >= 2:
            try:
                pts.append((float(row[0]), float(row[1])))
            except:
                pass
    return pts

def centre_gravite_polygone(data):
    """Centre de gravité d'un polygone."""
    pts = points_valides(data)
    n   = len(pts)
    if n < 3:
        return None, None
    A = Cx = Cy = 0.0
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i+1) % n]
        cross   = x1*y2 - x2*y1
        A  += cross
        Cx += (x1+x2)*cross
        Cy += (y1+y2)*cross
    A *= 0.5
    if A == 0:
        return None, None
    return Cx/(6*A), Cy/(6*A)

def aire_et_centre(data):
    """Aire et centre de gravité."""
    pts = points_valides(data)
    n   = len(pts)
    if n < 3:
        return 0.0, None, None
    A = Cx = Cy = 0.0
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i+1) % n]
        cross   = x1*y2 - x2*y1
        A  += cross
        Cx += (x1+x2)*cross
        Cy += (y1+y2)*cross
    A *= 0.5
    if A == 0:
        return 0.0, None, None
    return abs(A), Cx/(6*A), Cy/(6*A)

def aire_centre_inertie_origine(data):
    """Caractéristiques mécaniques / origine."""
    pts = points_valides(data)
    n   = len(pts)
    if n < 3:
        return 0.0, None, None, 0.0, 0.0
    A = Cx = Cy = Ix = Iy = 0.0
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i+1) % n]
        cross = x1*y2 - x2*y1
        A  += cross
        Cx += (x1+x2)*cross
        Cy += (y1+y2)*cross
        Ix += (y1**2 + y1*y2 + y2**2)*cross
        Iy += (x1**2 + x1*x2 + x2**2)*cross
    A *= 0.5
    if A == 0:
        return 0.0, None, None, 0.0, 0.0
    return np.abs(A), Cx/(6*A), Cy/(6*A), np.abs(Ix/12), np.abs(Iy/12)

def caracteristiques_mecaniques(contour, evidement):
    """Caractéristiques mécaniques avec évidements."""
    A, Cx, Cy, Ix0, Iy0 = aire_centre_inertie_origine(contour)
    if A == 0:
        return 0.0, None, None, 0.0, 0.0
    A_tot  = A
    Cx_tot = A*Cx
    Cy_tot = A*Cy
    Ix_tot = Ix0
    Iy_tot = Iy0
    N = n_liste(evidement)
    for i in range(1, N+1):
        Li = liste_i(evidement, i)
        Ai, cxi, cyi, ixi, iyi = aire_centre_inertie_origine(Li)
        if Ai > 0:
            A_tot  -= Ai
            Cx_tot -= Ai*cxi
            Cy_tot -= Ai*cyi
            Ix_tot -= ixi
            Iy_tot -= iyi
    if A_tot == 0:
        return 0.0, None, None, 0.0, 0.0
    CxG = Cx_tot/A_tot
    CyG = Cy_tot/A_tot
    return np.abs(A_tot), CxG, CyG, np.abs(Ix_tot - A_tot*CyG**2), np.abs(Iy_tot - A_tot*CxG**2)

def translation_points(data, Cx, Cy):
    """Translate des points vers le repère centré en (Cx, Cy)."""
    pts_trans = []
    for row in data:
        if len(row) >= 2:
            try:
                pts_trans.append((float(row[0])-Cx, float(row[1])-Cy))
            except:
                pass
    return pts_trans

def transformation_repere_cg(contour, evidement):
    """Transforme contour et évidements dans le repère CG."""
    A, Cx, Cy, Ix, Iy = caracteristiques_mecaniques(contour, evidement)
    if Cx is None or Cy is None:
        return None, [], None, None
    contour_cg    = translation_points(contour, Cx, Cy)
    evidements_cg = []
    N = n_liste(evidement)
    if N == 0:
        return contour_cg, evidements_cg, Cx, Cy
    for i in range(1, N+1):
        Li = liste_i(evidement, i)
        evidements_cg.append(translation_points(Li, Cx, Cy))
    return contour_cg, evidements_cg, Cx, Cy

def nettoyer_polygone(points, tol=1e-9):
    """Supprime le doublon de fermeture du polygone."""
    if len(points) < 2:
        return points
    x1, y1 = points[0]
    x2, y2 = points[-1]
    if abs(x1-x2) < tol and abs(y1-y2) < tol:
        return points[:-1]
    return points

def S_com(polygon1, evi, eps0, alpha, beta):
    """Surface comprimée de la section."""
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return 0.0
    def filter_comprime(x, y):
        eps     = eps0 + alpha*x + beta*y
        abs_eps = np.abs(eps)
        return (eps + abs_eps) / (2*abs_eps + 1e-15)
    Ic = polygone_integrate(filter_comprime, contour_cg)
    Iv = 0.0
    for trou in evidements_cg:
        trou = nettoyer_polygone(trou)
        if len(trou) >= 3:
            Iv += polygone_integrate(filter_comprime, trou)
    return Ic - Iv


# ══════════════════════════════════════════════════════════════════════
# ELS — INTÉGRATION APPROCHÉE (loi linéaire, pas de primitive simple)
# ══════════════════════════════════════════════════════════════════════

def calculer_N_My_Mz(polygon1, evi, p_acier, s_acier, a_com, n, eps0, alpha, beta):
    """Calcul N, My, Mz ELS — loi béton linéaire."""
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

    def f_complet(x, y):
        eps = eps0 + alpha*x + beta*y
        sig = sigma_c_n1(eps, n)
        return np.array([sig, sig*y, sig*x])

    res_c = polygone_integrate(f_complet, contour_cg)
    res_v = np.zeros(3)
    for trou in evidements_cg:
        if len(trou) >= 3:
            res_v += polygone_integrate(f_complet, trou)
    Nc, Mc_y, Mc_z = res_c - res_v

    Ns = Msy = Msz = 0.0
    if p_acier is not None and len(p_acier) > 0:
        acier = np.array(acier_G(p_acier, s_acier), dtype=float)
        if acier.size > 0:
            x_s    = acier[:, 0] - Cx
            y_s    = acier[:, 1] - Cy
            eps_s  = eps0 + alpha*x_s + beta*y_s
            sig_s  = sigma_s_lin1(eps_s, a_com)
            forces = sig_s * acier[:, 2]
            Ns     = np.sum(forces)
            Msy    = np.sum(forces * y_s)
            Msz    = np.sum(forces * x_s)

    return Ns+Nc, Msy+Mc_y, Msz+Mc_z

def solve_GG_ELS(polygon1, evi, p_acier, s_acier, a_com, n, Nobj, Myobj, Mzobj):
    """Résolution ELS."""
    x0 = np.array([0., 0.001, 0.0])
    def residuals(eps):
        N, My, Mz = calculer_N_My_Mz(polygon1, evi, p_acier, s_acier, a_com, n, *eps)
        return np.array([N-Nobj, My-Myobj, Mz-Mzobj])
    def jacobian(eps, h=1e-6):
        J  = np.empty((3, 3))
        r0 = residuals(eps)
        for i in range(3):
            dh = np.zeros(3); dh[i] = h
            J[:, i] = (residuals(eps+dh) - r0) / h
        return J
    result = root(residuals, x0, jac=jacobian, method="hybr",
                  tol=1e-5, options={"maxfev": 800})
    return result.x

def resultats_GG_ELS(polygon1, evi, p_acier, s_acier, a_com, n, eps0, alpha, beta):
    """Résultats complets ELS."""
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return {}
    poly_array = np.asarray(contour_cg, dtype=float)
    eps_c     = eps0 + alpha*poly_array[:, 0] + beta*poly_array[:, 1]
    espcmax   = float(np.max(eps_c))
    espcmin   = float(np.min(eps_c))
    sig_c_max = float(sigma_c_n1(espcmax, n))
    sig_c_min = float(sigma_c_n1(espcmin, n))

    espsmax = espsmin = sig_s_max = sig_s_min = ""
    if p_acier is not None and len(p_acier) > 0:
        acier = acier_G(p_acier, s_acier)
        if len(acier) > 0:
            acier_array = np.asarray(acier, dtype=float)
            x_s   = acier_array[:, 0] - Cx
            y_s   = acier_array[:, 1] - Cy
            eps_s = eps0 + alpha*x_s + beta*y_s
            espsmax   = float(np.max(eps_s))
            espsmin   = float(np.min(eps_s))
            sig_s_max = float(sigma_s_lin1(espsmax, a_com))
            sig_s_min = float(sigma_s_lin1(espsmin, a_com))

    A_com    = S_com(polygon1, evi, eps0, alpha, beta)
    N, My, Mz = calculer_N_My_Mz(polygon1, evi, p_acier, s_acier, a_com, n, eps0, alpha, beta)

    return {
        "ACOM": A_com,
        "EPS_C_MAX": espcmax,  "EPS_C_MIN": espcmin,
        "SIG_C_MAX": sig_c_max, "SIG_C_MIN": sig_c_min,
        "EPS_S_MAX": espsmax,  "EPS_S_MIN": espsmin,
        "SIG_S_MAX": sig_s_max, "SIG_S_MIN": sig_s_min,
        "N": N, "MY": My, "MZ": Mz, "CX": Cx, "CY": Cy
    }

def e_resultats_GG_ELS(polygon1, evi, p_acier, s_acier, a_com, n, eps0, alpha, beta, resultats):
    """Extraction sélective résultats ELS."""
    tout = resultats_GG_ELS(polygon1, evi, p_acier, s_acier, a_com, n, eps0, alpha, beta)
    return [tout[r] for r in resultats.split(',')]


# ══════════════════════════════════════════════════════════════════════
# ELU — INTÉGRATION ANALYTIQUE EXACTE PAR GREEN (n=2, C20-C50)
# ══════════════════════════════════════════════════════════════════════

def calculer_N_My_Mz_ELU_pararect(
    polygon1, evi, p_acier, s_acier,
    a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
    eps0, alpha, beta
):
    """
    Calcul ELU N, My, Mz — intégration ANALYTIQUE EXACTE par Green.
    Valide pour béton C20-C50 (n=2).
    Indépendant du sens de saisie du polygone.
    """
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return 0.0, 0.0, 0.0

    e2 = float(eps_c2(fck))

    # ── Béton plein ─────────────────────────────────────────────────
    Nc, Mc_y, Mc_z = _integrer_polygone_exact(
        contour_cg, fcd, e2, eps0, alpha, beta
    )

    # ── Évidements ──────────────────────────────────────────────────
    for trou in evidements_cg:
        trou = nettoyer_polygone(trou)
        if len(trou) >= 3:
            dN, dMy, dMz = _integrer_polygone_exact(
                trou, fcd, e2, eps0, alpha, beta
            )
            Nc   -= dN
            Mc_y -= dMy
            Mc_z -= dMz

    # ── Acier ────────────────────────────────────────────────────────
    Ns = Msy = Msz = 0.0
    if p_acier is not None and len(p_acier) > 0:
        acier = np.array(acier_G(p_acier, s_acier), dtype=float)
        if acier.size > 0:
            x_s    = acier[:, 0] - Cx
            y_s    = acier[:, 1] - Cy
            eps_s  = eps0 + alpha*x_s + beta*y_s
            sig_s  = sigma_s_palier1(fyd, k, eps_uk, eps_ud, eps_s, a_com)
            forces = sig_s * acier[:, 2]
            Ns     = np.sum(forces)
            Msy    = np.sum(forces * y_s)
            Msz    = np.sum(forces * x_s)

    return Ns+Nc, Msy+Mc_y, Msz+Mc_z


def solve_GG_ELU_pararect(
    polygon1, evi, p_acier, s_acier,
    a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
    Nobj, Myobj, Mzobj
):
    """Résolution ELU — intégration analytique exacte."""
    x0 = np.array([0.0, 0.1, 0.0])

    def residuals(eps):
        N, My, Mz = calculer_N_My_Mz_ELU_pararect(
            polygon1, evi, p_acier, s_acier,
            a_com, fck, fcd, fyd, k, eps_uk, eps_ud, *eps
        )
        return np.array([N-Nobj, My-Myobj, Mz-Mzobj])

    def jacobian(eps, h=1e-6):
        J  = np.empty((3, 3))
        r0 = residuals(eps)
        for i in range(3):
            dh = np.zeros(3); dh[i] = h
            J[:, i] = (residuals(eps+dh) - r0) / h
        return J

    result = root(residuals, x0, jac=jacobian, method="hybr",
                  tol=1e-6, options={"maxfev": 1000})
    return result.x


def resultats_GG_ELU_pararect(
    polygon1, evi, p_acier, s_acier,
    a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
    eps0, alpha, beta
):
    """Résultats complets ELU — intégration analytique exacte."""
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
    if contour_cg is None:
        return {}

    poly_array = np.asarray(contour_cg, dtype=float)
    eps_c      = eps0 + alpha*poly_array[:, 0] + beta*poly_array[:, 1]
    espcmax    = float(np.max(eps_c))
    espcmin    = float(np.min(eps_c))
    sig_c_max  = float(sigma_c_pararect1(fck, fcd, espcmax))
    sig_c_min  = float(sigma_c_pararect1(fck, fcd, espcmin))

    epssmax = epssmin = sig_s_max = sig_s_min = ""
    pour_a_pl = 0.0
    if p_acier is not None and len(p_acier) > 0:
        acier = acier_G(p_acier, s_acier)
        if len(acier) > 0:
            acier_array = np.asarray(acier, dtype=float)
            x_s   = acier_array[:, 0] - Cx
            y_s   = acier_array[:, 1] - Cy
            eps_s = eps0 + alpha*x_s + beta*y_s
            epssmax   = float(np.max(eps_s))
            epssmin   = float(np.min(eps_s))
            plast     = np.sum(np.abs(eps_s) >= (float(fyd) / 200.0))
            pour_a_pl = float(plast / len(eps_s) * 100.0)
            sig_s_max = float(sigma_s_palier1(fyd, k, eps_uk, eps_ud, epssmax, a_com))
            sig_s_min = float(sigma_s_palier1(fyd, k, eps_uk, eps_ud, epssmin, a_com))

    A_com     = float(S_com(polygon1, evi, eps0, alpha, beta))
    N, My, Mz = calculer_N_My_Mz_ELU_pararect(
        polygon1, evi, p_acier, s_acier,
        a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
        eps0, alpha, beta
    )

    return {
        "EPS0": eps0,    "ALPHA": alpha,  "BETA": beta,
        "ACOM": A_com,
        "EPS_C_MAX": espcmax,   "EPS_C_MIN": espcmin,
        "SIG_C_MAX": sig_c_max, "SIG_C_MIN": sig_c_min,
        "EPS_S_MAX": epssmax,   "EPS_S_MIN": epssmin,
        "SIG_S_MAX": sig_s_max, "SIG_S_MIN": sig_s_min,
        "N": float(N), "MY": float(My), "MZ": float(Mz),
        "PA": pour_a_pl, "CX": Cx, "CY": Cy
    }


def e_resultats_GG_ELU_pararect(
    polygon1, evi, p_acier, s_acier,
    a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
    eps0, alpha, beta, resultats
):
    """Extraction sélective résultats ELU."""
    tout = resultats_GG_ELU_pararect(
        polygon1, evi, p_acier, s_acier,
        a_com, fck, fcd, fyd, k, eps_uk, eps_ud,
        eps0, alpha, beta
    )
    return [tout[r] for r in resultats.split(',')]
