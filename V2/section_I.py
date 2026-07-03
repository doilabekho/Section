import numpy as np
import scipy.integrate 
from scipy.integrate import quad
from scipy.optimize import root, fsolve, least_squares, root_scalar
from scipy.spatial import Delaunay
from matplotlib.path import Path
from xlwings import func
from scipy.special import roots_legendre

from beton import *   # ← fonctionne en local et sur GitHub
from acier import *   # ← fonctionne en local et sur GitHub	

_ES =200000
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
# ════════════════════════════════════════════════════════════════════════════
# 3. GÉOMÉTRIE — SECTION EN I
# ════════════════════════════════════════════════════════════════════════════
#
#  Paramètres géométriques (tous en mm) :
#    b, h      : largeur et hauteur totale de la section
#    bs, hs    : largeur et hauteur de la table supérieure
#    gs        : congé supérieur (gousset)
#    bi, hi    : largeur et hauteur de la table inférieure
#    gi        : congé inférieur (gousset)
#
#   ┌──────────────── bs ─────────────────┐
#   │                                     │  hs
#   └──────┬─────────────────────┬────────┘
#          │ gs                  │ gs
#          │◄── b/2              │
#          │                     │
#          │       âme           │
#          │                     │
#          │ gi                  │ gi
#   ┌──────┴─────────────────────┴────────┐
#   │                                     │  hi
#   └──────────────── bi ─────────────────┘

# ════════════════════════════════════════════════════════════════════════════
# 3. GÉOMÉTRIE — SECTION EN I
# ════════════════════════════════════════════════════════════════════════════
#
#  Paramètres géométriques (tous en mm) :
#    b, h      : largeur et hauteur totale de la section
#    bs, hs    : largeur et hauteur de la table supérieure
#    gs        : congé supérieur (gousset)
#    bi, hi    : largeur et hauteur de la table inférieure
#    gi        : congé inférieur (gousset)
#
#   ┌──────────────── bs ─────────────────┐
#   │                                     │  hs
#   └──────┬─────────────────────┬────────┘
#          │ gs                  │ gs
#          │◄── b/2              │
#          │                     │
#          │       âme           │
#          │                     │
#          │ gi                  │ gi
#   ┌──────┴─────────────────────┴────────┐
#   │                                     │  hi
#   └──────────────── bi ─────────────────┘

# ════════════════════════════════════════════════════════════════════════════
# 1. GÉOMÉTRIE ET POLYNÔMES (Repère initial et centré)
# ════════════════════════════════════════════════════════════════════════════

def trans_i(b, h, bs, hs, gs, bi, hi, gi):
    """Retourne les sommets du polygone section I (sens trigonométrique)."""
    return np.array([
        [-bi/2, 0],       [-bi/2+bi, 0],
        [ bi/2, hi],      [ b/2,     hi + gi],
        [ b/2,  h-hs-gs], [ bs/2,    h - hs],
        [ bs/2, h],       [-bs/2,    h],
        [-bs/2, h-hs],    [-b/2,     h-hs-gs],
        [-b/2,  hi+gi],   [-bi/2,    hi],
        [-bi/2, 0],
    ], dtype=float)


def Nc_Gy_ELS(b, h, bs, hs, gs, bi, hi, gi):
    """Ordonnée du centre de gravité (fibre inf = 0)."""
    pts = trans_i(b, h, bs, hs, gs, bi, hi, gi)
    x0, y0 = pts[:-1, 0], pts[:-1, 1]
    x1, y1 = pts[1:,  0], pts[1:,  1]
    aire_el = x0 * y1 - x1 * y0
    aire    = 0.5 * aire_el.sum()
    yG      = (((y0 + y1) * aire_el).sum()) / (6.0 * aire)
    return float(yG)


def section_I(b, h, bs, hs, gs, bi, hi, gi):
    """Polygone centré sur le CDG (y_CDG = 0)."""
    yG   = Nc_Gy_ELS(b, h, bs, hs, gs, bi, hi, gi)
    poly = trans_i(b, h, bs, hs, gs, bi, hi, gi)
    poly[:, 1] -= yG
    return poly


# ════════════════════════════════════════════════════════════════════════════
# 2. MOTEUR D'INTÉGRATION UNIFIÉ (Centré analytiquement sur le CDG)
# ════════════════════════════════════════════════════════════════════════════

import numpy as np

def calculer_NM_beton(pts, eps0, beta, mode='ELS', n_els=15.0, fck=30.0, fcd=20.0):
    """
    Calcule [Nc, Mc_y] en une seule passe à partir du polygone 'pts' centré au CDG.
    
    Parameters:
    -----------
    pts   : np.ndarray des sommets de la section (généré par section_I)
    eps0  : déformation au CDG (y = 0)
    beta  : courbure (gradient de déformation selon y)
    mode  : 'ELS' ou 'ELU'
    """
    # 1. Extraction directe des altitudes (y) des 5 couches depuis 'pts'
    # 'pts' étant centré sur le CDG, ces valeurs intègrent déjà le décalage -yG
    y0 = pts[0, 1]   # Base de la table inférieure
    y1 = pts[2, 1]   # Sommet de la table inférieure / Bas du gousset inf
    y2 = pts[3, 1]   # Sommet du gousset inf / Bas de l'âme
    y3 = pts[4, 1]   # Sommet de l'âme / Bas du gousset sup
    y4 = pts[5, 1]   # Sommet du gousset sup / Bas de la table sup
    y5 = pts[6, 1]   # Sommet de la table supérieure
    
    # 2. Extraction des largeurs (b) par différence des coordonnées X (X_droit - X_gauche)
    bi_start = pts[1, 0] - pts[0, 0]  # Largeur bi de la table inf
    bi_end   = pts[2, 0] - pts[11, 0] # Fin de la table inf
    b_start  = pts[3, 0] - pts[10, 0] # Largeur b de l'âme
    b_end    = pts[4, 0] - pts[9, 0]  # Fin de l'âme
    bs_start = pts[5, 0] - pts[8, 0]  # Largeur bs de la table sup
    bs_end   = pts[6, 0] - pts[7, 0]  # Sommet de la table sup
    
    # Définition des 5 couches géométriques : (y_min, y_max, b_départ, b_arrivée)
    couches = [
        (y0, y1, bi_start, bi_end),  # 1. Table inférieure (Rectangle)
        (y1, y2, bi_end, b_start),   # 2. Gousset inférieur (Trapèze)
        (y2, y3, b_start, b_end),    # 3. Âme (Rectangle)
        (y3, y4, b_end, bs_start),   # 4. Gousset supérieur (Trapèze)
        (y4, y5, bs_start, bs_end)   # 5. Table supérieure (Rectangle)
    ]
    
    # 3. Paramétrage des constantes matériaux selon le mode
    if mode == 'ELS':
        C_els = 200000.0 / (float(n_els) * 1000.0)
        e2, expo_elu = 0.0, 2.0
        n_gauss = 3
    else:  # mode == 'ELU'
        e2 = eps_c2(fck)
        expo_elu = eps_n(fck)
        n_gauss = 5

    gauss_points, gauss_weights = roots_legendre(n_gauss)

    
    Nc = 0.0
    Mc_y = 0.0
    
    # 4. Boucle d'intégration
    for y_start, y_end, b_start, b_end in couches:
        if y_start >= y_end:
            continue
            
        splits = [y_start, y_end]
        if abs(beta) > 1e-12:
            y_zero = -eps0 / beta
            if y_start < y_zero < y_end:
                splits.append(y_zero)
                
            if mode == 'ELU':
                y_e2 = (e2 - eps0) / beta
                if y_start < y_e2 < y_end:
                    splits.append(y_e2)
                    
        splits = sorted(list(set(splits)))
        
        for i in range(len(splits) - 1):
            y_min, y_max = splits[i], splits[i+1]
            half_len = (y_max - y_min) / 2.0
            mean_val = (y_max + y_min) / 2.0
            
            for xi, wi in zip(gauss_points, gauss_weights):
                y_g = half_len * xi + mean_val
                eps_g = eps0 + beta * y_g
                
                # Calcul de la contrainte sigma
                if eps_g <= 0:
                    sig = 0.0
                elif mode == 'ELS':
                    sig = C_els * eps_g
                else: 
                    sig = fcd * (1.0 - (1.0 - eps_g / e2)**expo_elu) if eps_g <= e2 else fcd
                
                # Interpolation de la largeur b(y) au point de Gauss
                t = (y_g - y_start) / (y_end - y_start) if abs(y_end - y_start) > 1e-11 else 0.0
                b_g = b_start + t * (b_end - b_start)
                
                facteur_integration = wi * half_len * b_g
                
                Nc += sig * facteur_integration
                Mc_y += sig * y_g * facteur_integration
                
    return Nc, Mc_y

# ════════════════════════════════════════════════════════════════════════════
# 3. INTERFACES COMPORTEMENTALES (Passage obligé par le dictionnaire 'geom')
# ════════════════════════════════════════════════════════════════════════════

def _pts_I(b, h, bs, hs, gs, bi, hi, gi):
    """Raccourci pour obtenir le polygone centré."""
    return section_I(b, h, bs, hs, gs, bi, hi, gi)


def _yG(b, h, bs, hs, gs, bi, hi, gi):
    """Raccourci pour obtenir l'ordonnée absolue du CDG."""
    return Nc_Gy_ELS(b, h, bs, hs, gs, bi, hi, gi)



def _NM_beton_ELS(pts, n, eps0, beta):
    """
    Calcule [Nc, Mc] en ELS (UNE seule passe via Gauss 1D).
    Le premier argument est désormais le dictionnaire 'geom'.
    """
    return calculer_NM_beton(pts, eps0, beta, mode='ELS', n_els=n)




def _NM_acier_ELS(yG, h, asup, ainf, esup, einf, eps0, beta):
    """Efforts acier — ELS loi linéaire."""
    ys_sup =  h - yG - esup   # ordonnée acier sup / CDG
    ys_inf = -yG      + einf  # ordonnée acier inf / CDG
    sig_sup = sigma_s_lin(eps0 + beta * ys_sup)
    sig_inf = sigma_s_lin(eps0 + beta * ys_inf)
    Ns  = sig_sup * asup * 1e-4 + sig_inf * ainf * 1e-4
    Ms  = sig_sup * asup * 1e-4 * ys_sup + sig_inf * ainf * 1e-4 * ys_inf
    return Ns, Ms


@func
def S_com(b, h, bs, hs, gs, bi, hi, gi, eps0, beta):
    """
    Surface comprimée (zone où ε > 0).
    Calcul analytique exact par couches, sans polygone_integrate ni approximation.
    """
    # 1. Génération du polygone centré au CDG
    pts = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    
    # 2. Extraction des altitudes (y) des 5 couches (déjà centrées sur le CDG)
    y0 = pts[0, 1]   # Base de la table inférieure
    y1 = pts[2, 1]   # Sommet de la table inférieure / Bas du gousset inf
    y2 = pts[3, 1]   # Sommet du gousset inf / Bas de l'âme
    y3 = pts[4, 1]   # Sommet de l'âme / Bas du gousset sup
    y4 = pts[5, 1]   # Sommet du gousset sup / Bas de la table sup
    y5 = pts[6, 1]   # Sommet de la table supérieure
    
    # 3. Extraction des largeurs (b) associées
    bi_start = pts[1, 0] - pts[0, 0]  # Largeur bi de la table inf
    bi_end   = pts[2, 0] - pts[11, 0] # Fin de la table inf
    b_start  = pts[3, 0] - pts[10, 0] # Largeur b de l'âme
    b_end    = pts[4, 0] - pts[9, 0]  # Fin de l'âme
    bs_start = pts[5, 0] - pts[8, 0]  # Largeur bs de la table sup
    bs_end   = pts[6, 0] - pts[7, 0]  # Sommet de la table sup
    
    couches = [
        (y0, y1, bi_start, bi_end),  # 1. Table inférieure
        (y1, y2, bi_end, b_start),   # 2. Gousset inférieur
        (y2, y3, b_start, b_end),    # 3. Âme
        (y3, y4, b_end, bs_start),   # 4. Gousset supérieur
        (y4, y5, bs_start, bs_end)   # 5. Table supérieure
    ]
    
    s_com_total = 0.0
    
    # 4. Boucle de détection et de sommation des aires comprimées
    for y_start, y_end, b_start, b_end in couches:
        if y_start >= y_end:
            continue
            
        # Trouver la frontière de fissuration (ε = 0) au sein de la couche
        splits = [y_start, y_end]
        if abs(beta) > 1e-12:
            y_zero = -eps0 / beta
            if y_start < y_zero < y_end:
                splits.append(y_zero)
                
        splits = sorted(list(set(splits)))
        
        # Intégration géométrique des sous-intervalles
        for i in range(len(splits) - 1):
            y_min, y_max = splits[i], splits[i+1]
            
            # Test du signe de la déformation au milieu du sous-intervalle
            y_milieu = (y_min + y_max) / 2.0
            eps_milieu = eps0 + beta * y_milieu
            
            # Si le segment est comprimé (ε > 0), on ajoute son aire géométrique
            if eps_milieu > 0:
                # Interpolation locale des largeurs aux bornes y_min et y_max
                if abs(y_end - y_start) > 1e-11:
                    t_min = (y_min - y_start) / (y_end - y_start)
                    t_max = (y_max - y_start) / (y_end - y_start)
                    b_min = b_start + t_min * (b_end - b_start)
                    b_max = b_start + t_max * (b_end - b_start)
                else:
                    b_min, b_max = b_start, b_end
                
                # Formule analytique exacte de l'aire d'un trapèze : (b1 + b2) * h / 2
                aire_trapeze = 0.5 * (b_min + b_max) * (y_max - y_min)
                s_com_total += aire_trapeze
                
    return float(s_com_total)


@func
def Nc_I_ELS(b, h, bs, hs, gs, bi, hi, gi, n, eps0, beta):
    """Effort normal béton — ELS."""
    pts = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    Nc, _ = _NM_beton_ELS(pts, n, eps0, beta)
    return float(Nc)


@func
def Mc_I_ELS(b, h, bs, hs, gs, bi, hi, gi, n, eps0, beta):
    """Moment béton / CDG — ELS."""
    pts = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    _, Mc = _NM_beton_ELS(pts, n, eps0, beta)
    return float(Mc)


@func
def Ns_I_ELS(b, h, bs, hs, gs, bi, hi, gi, asup, ainf, esup, einf, n, eps0, beta):
    """Effort normal acier — ELS."""
    yG       = _yG(b, h, bs, hs, gs, bi, hi, gi)
    Ns, _    = _NM_acier_ELS(yG, h, asup, ainf, esup, einf, eps0, beta)
    return float(Ns)


@func
def Ms_I_ELS(b, h, bs, hs, gs, bi, hi, gi, asup, ainf, esup, einf, n, eps0, beta):
    """Moment acier / CDG — ELS."""
    yG       = _yG(b, h, bs, hs, gs, bi, hi, gi)
    _, Ms    = _NM_acier_ELS(yG, h, asup, ainf, esup, einf, eps0, beta)
    return float(Ms)


@func
def N_I_ELS(b, h, bs, hs, gs, bi, hi, gi, asup, ainf, esup, einf, n, eps0, beta):
    """Effort normal total — ELS."""
    pts      = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    yG       = _yG(b, h, bs, hs, gs, bi, hi, gi)
    Nc, _    = _NM_beton_ELS(pts, n, eps0, beta)
    Ns, _    = _NM_acier_ELS(yG, h, asup, ainf, esup, einf, eps0, beta)
    return float(Nc + Ns)


@func
def M_I_ELS(b, h, bs, hs, gs, bi, hi, gi, asup, ainf, esup, einf, n, eps0, beta):
    """Moment total / CDG — ELS."""
    pts      = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    yG       = _yG(b, h, bs, hs, gs, bi, hi, gi)
    _, Mc    = _NM_beton_ELS(pts, n, eps0, beta)
    _, Ms    = _NM_acier_ELS(yG, h, asup, ainf, esup, einf, eps0, beta)
    return float(Mc + Ms)


# ════════════════════════════════════════════════════════════════════════════
# SOLVEUR ELS — (N, M) → (ε₀, β)
# ════════════════════════════════════════════════════════════════════════════

def _residuals_ELS(ep, b, h, bs, hs, gs, bi, hi, gi,
                    asup, ainf, esup, einf, n,
                    Nobj, Mobj):
    """Résidu [N−Nobj, M−Mobj] pour un vecteur ep = [ε₀, β]."""
    N, M = calculer_N_M(
        b, h, bs, hs, gs, bi, hi, gi,
        asup, ainf, esup, einf, n,
        ep[0], ep[1]
    )
    return np.array([N - Nobj, M - Mobj])


def _grille_x0_ELS(Nobj, Mobj, h):
    """
    Génère des points de départ couvrant le domaine ELS :
      ε₀ ∈ [−2.5, +2.5] ‰   (traction pure → compression modérée, acier limite ~fyk/Es)
      β  ∈ [−3.5/h, +3.5/h] ‰/mm
    """
    eps0_vals = [-2.0, -0.5, 0.0, 0.5, 1.5]
    beta_vals = [-2.5/h, 0.0, 2.5/h]
    beta_est  = 2.5/h if Mobj > 0 else -2.5/h
    eps0_est  = 0.0
    points = [(eps0_est, beta_est)]
    for e in eps0_vals:
        for b_ in beta_vals:
            points.append((e, b_))
    return [np.array([e, b_]) for e, b_ in points]


def _solve_hybr_ELS(resid_fn, x0):
    try:
        sol = root(
            resid_fn, x0,
            jac=lambda ep: _jacobian_centre(ep, resid_fn),
            method="hybr",
            tol=1e-6,
            options={"maxfev": 600}
        )
        if sol.success and np.max(np.abs(sol.fun)) < 1e-3:
            return sol.x
    except Exception:
        pass
    return None


def _solve_trf_ELS(resid_fn, x0):
    try:
        sol = least_squares(
            resid_fn, x0,
            jac=lambda ep: _jacobian_centre(ep, resid_fn),
            method='trf',
            xtol=1e-6, ftol=1e-6, gtol=1e-6,
            max_nfev=400
        )
        if np.max(np.abs(sol.fun)) < 1e-3:
            return sol.x
    except Exception:
        pass
    return None


def _solve_multistart_ELS(resid_fn, Nobj, Mobj, h):
    best_x, best_err = None, np.inf
    for x0 in _grille_x0_ELS(Nobj, Mobj, h):
        x = _solve_hybr_ELS(resid_fn, x0)
        if x is None:
            x = _solve_trf_ELS(resid_fn, x0)
        if x is not None:
            err = float(np.max(np.abs(resid_fn(x))))
            if err < best_err:
                best_err, best_x = err, x
    if best_err < 1e-2:
        return best_x
    return None


def _solve_bissection_ELS(resid_fn, Nobj, Mobj, h):
    """Filet de sécurité ultime : bissection séquentielle β puis ε₀."""
    from scipy.optimize import brentq as _brentq
    eps_lim = 2.5  # ‰ — plage large pour rester couvrant

    def beta_from_eps0(eps0):
        def fM(beta):
            return float(resid_fn(np.array([eps0, beta]))[1])
        try:
            return _brentq(fM, -eps_lim/h, eps_lim/h, xtol=1e-8, maxiter=100)
        except ValueError:
            return None

    def fN(eps0):
        beta = beta_from_eps0(eps0)
        if beta is None:
            return np.nan
        return float(resid_fn(np.array([eps0, beta]))[0])

    eps0_grid = np.linspace(-eps_lim, eps_lim, 30)
    fN_vals = [fN(e) for e in eps0_grid]

    for i in range(len(fN_vals) - 1):
        f1, f2 = fN_vals[i], fN_vals[i+1]
        if np.isnan(f1) or np.isnan(f2):
            continue
        if f1 * f2 < 0:
            try:
                eps0_sol = _brentq(fN, eps0_grid[i], eps0_grid[i+1],
                                    xtol=1e-7, maxiter=80)
                beta_sol = beta_from_eps0(eps0_sol)
                if beta_sol is not None:
                    x = np.array([eps0_sol, beta_sol])
                    if np.max(np.abs(resid_fn(x))) < 1e-2:
                        return x
            except Exception:
                continue
    return None


@func
def solve_I_ELS(b, h, bs, hs, gs, bi, hi, gi, asup, ainf, esup, einf, n, Nobj, Mobj):
    """
    Résout (N = Nobj, M = Mobj) → (ε₀, β) en ELS.

    Ordre de résolution :
      0. Cas analytique "section entièrement tendue" (fermé, instantané)
      1. root/hybr
      2. least_squares/trf
      3. Multi-start (grille de points de départ)
      4. Bissection 1D (filet ultime)

    Retourne [ε₀, β] ou [nan, nan] si aucune méthode ne converge.
    Convention : eps(y) = eps0 + beta*(y - yG), y positif vers le bas depuis la fibre sup.
    """

    yG = Nc_Gy_ELS(b, h, bs, hs, gs, bi, hi, gi)

    cond_M_pos = (Mobj > 0 and Mobj / (Nobj - 1e-15) <= -(yG - einf - esup))

    # ── Étape 0 : cas analytique (section tout tendue) ─────────────────────
    if Nobj < 0 and asup * ainf > 0 and cond_M_pos:
        y_sup = esup
        y_inf = h - einf
        a = y_sup - yG
        c = y_inf - yG
        d = a - c
        if abs(d) > 1e-9:
            F_sup = (Mobj - Nobj * c) / d
            F_inf = Nobj - F_sup

            sigma_sup = F_sup / asup * 10000
            sigma_inf = F_inf / ainf * 10000

            eps_sup = sigma_sup / 200000
            eps_inf = sigma_inf / 200000

            beta = (eps_sup - eps_inf) / (a - c)
            eps0 = eps_sup - beta * a

            sol = np.array([eps0, beta])
            N_chk, M_chk = calculer_N_M(b, h, bs, hs, gs, bi, hi, gi,
                                         asup, ainf, esup, einf, n, *sol)
            if (abs(N_chk - Nobj) < 1e-3 * max(1, abs(Nobj)) and
                    abs(M_chk - Mobj) < 1e-3 * max(1, abs(Mobj))):
                return sol
            # sinon on continue vers la cascade numérique ci-dessous

    # ── Cascade numérique ────────────────────────────────────────────────
    def resid(ep):
        return _residuals_ELS(
            ep, b, h, bs, hs, gs, bi, hi, gi,
            asup, ainf, esup, einf, n,
            Nobj, Mobj
        )

    x0 = np.array([0., 0.001])

    x = _solve_hybr_ELS(resid, x0)
    if x is not None:
        return x

    x = _solve_trf_ELS(resid, x0)
    if x is not None:
        return x

    x = _solve_multistart_ELS(resid, Nobj, Mobj, h)
    if x is not None:
        return x

    x = _solve_bissection_ELS(resid, Nobj, Mobj, h)
    if x is not None:
        return x

    return np.array([np.nan, np.nan])

# ════════════════════════════════════════════════════════════════════════════
# 7. Résultats ELS — Calculs contrainte
# ════════════════════════════════════════════════════════════════════════════    
@func
def resultats_I_ELS(b, h, bs, hs, gs, bi, hi, gi, asup, ainf, esup, einf, n, eps0, beta):
    yc = Nc_Gy_ELS(b, h, bs, hs, gs, bi, hi, gi)
    epscsup = eps0 + (h - yc) * beta
    epscinf = eps0 - yc * beta
    epsssup = eps0 + (h - yc - esup) * beta
    epssinf = eps0 - (yc - einf) * beta
    sigmacsup = sigma_c_n(epscsup, n)
    sigmacinf = sigma_c_n(epscinf, n)
    sigmassup = sigma_s_lin(epsssup)
    sigmasinf = sigma_s_lin(epssinf)

    N, M = calculer_N_M(b, h, bs, hs, gs, bi, hi, gi, asup, ainf, esup, einf, n, eps0, beta)
    # Calcul de la hauteur comprimée
   
    ycomp = lambda y: np.where((eps0 + beta * y) > 0, 1.0, 0.0)
    hcomp = quad(ycomp, -yc, h - yc)[0]

    if hcomp == h:
        etat = "entièrement comprimé"
    elif hcomp > 0:
        etat = "partiellement tendu"
    else:
        etat = "entièrement tendu"

    res = {
        "EPS0": eps0,
        "BETA": beta,
        "H_compr": hcomp,
        "Etat": etat,
        "EPS_C_SUP": epscsup,
        "EPS_S_SUP": epsssup,
        "EPS_S_INF": epssinf,
        "EPS_C_INF": epscinf,
        "SIG_C_SUP": sigmacsup,
        "SIG_S_SUP": sigmassup,
        "SIG_S_INF": sigmasinf,
        "SIG_C_INF": sigmacinf,
        "N": N,
        "M": M
    }
    return res
@func
def e_resultats_I_ELS(b, h, bs, hs, gs, bi, hi, gi, asup, ainf, esup, einf, n, eps0, beta,resultats):
    resultats_Excel = []
    resultats_tout = resultats_I_ELS(b, h, bs, hs, gs, bi, hi, gi, asup, ainf, esup, einf, n, eps0, beta)                    
    resultats_list = resultats.split(',')
    for r in resultats_list:       
        resultats_Excel.append(resultats_tout[r])
    return resultats_Excel   # résultats en ligne

# ════════════════════════════════════════════════════════════════════════════
# 8. EFFORTS INTERNES — ELU  (loi parabole-rectangle)
# ════════════════════════════════════════════════════════════════════════════

def _NM_beton_ELU(pts, fck, fcd, eps0, beta):
    """
    Calcule [Nc, Mc] en ELU (UNE seule passe via Gauss 1D).
    """
    return calculer_NM_beton(pts, eps0, beta, mode='ELU', fck=fck, fcd=fcd)


def _NM_acier_ELU(yG, h, asup, ainf, esup, einf,
                  fyd, k, eps_uk, eps_ud, eps0, beta):
    """Efforts acier — ELU loi palier."""
    ys_sup =  h - yG - esup
    ys_inf = -yG      + einf
    sig_sup = sigma_s_palier(fyd, k, eps_uk, eps_ud, eps0 + beta * ys_sup)
    sig_inf = sigma_s_palier(fyd, k, eps_uk, eps_ud, eps0 + beta * ys_inf)
    Ns  = sig_sup * asup * 1e-4 + sig_inf * ainf * 1e-4
    Ms  = sig_sup * asup * 1e-4 * ys_sup + sig_inf * ainf * 1e-4 * ys_inf
    return float(Ns), float(Ms)


@func
def Nc_I_ELU_pararect(b, h, bs, hs, gs, bi, hi, gi, fck, fcd, eps0, beta):
    """Effort normal béton — ELU."""
    pts = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    Nc, _ = _NM_beton_ELU(pts, fck, fcd, eps0, beta)
    return float(Nc)


@func
def Mc_I_ELU_pararect(b, h, bs, hs, gs, bi, hi, gi, fck, fcd, eps0, beta):
    """Moment béton / CDG — ELU."""
    pts = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    _, Mc = _NM_beton_ELU(pts, fck, fcd, eps0, beta)
    return float(Mc)


@func
def Ns_I_ELU_pararect(b, h, bs, hs, gs, bi, hi, gi,
                       asup, ainf, esup, einf,
                       fyd, k, eps_uk, eps_ud, eps0, beta):
    """Effort normal acier — ELU."""
    yG    = _yG(b, h, bs, hs, gs, bi, hi, gi)
    Ns, _ = _NM_acier_ELU(yG, h, asup, ainf, esup, einf,
                           fyd, k, eps_uk, eps_ud, eps0, beta)
    return Ns


@func
def Ms_I_ELU_pararect(b, h, bs, hs, gs, bi, hi, gi,
                       asup, ainf, esup, einf,
                       fyd, k, eps_uk, eps_ud, eps0, beta):
    """Moment acier / CDG — ELU."""
    yG    = _yG(b, h, bs, hs, gs, bi, hi, gi)
    _, Ms = _NM_acier_ELU(yG, h, asup, ainf, esup, einf,
                           fyd, k, eps_uk, eps_ud, eps0, beta)
    return Ms


@func
def N_I_ELU_pararect(b, h, bs, hs, gs, bi, hi, gi,
                      asup, ainf, esup, einf,
                      fck, fcd, fyd, k, eps_uk, eps_ud, eps0, beta):
    """Effort normal total — ELU."""
    pts      = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    yG       = _yG(b, h, bs, hs, gs, bi, hi, gi)
    Nc, _    = _NM_beton_ELU(pts, fck, fcd, eps0, beta)
    Ns, _    = _NM_acier_ELU(yG, h, asup, ainf, esup, einf,
                              fyd, k, eps_uk, eps_ud, eps0, beta)
    return float(Nc + Ns)


@func
def M_I_ELU_pararect(b, h, bs, hs, gs, bi, hi, gi,
                      asup, ainf, esup, einf,
                      fck, fcd, fyd, k, eps_uk, eps_ud, eps0, beta):
    """Moment total / CDG — ELU."""
    pts      = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    yG       = _yG(b, h, bs, hs, gs, bi, hi, gi)
    _, Mc    = _NM_beton_ELU(pts, fck, fcd, eps0, beta)
    _, Ms    = _NM_acier_ELU(yG, h, asup, ainf, esup, einf,
                              fyd, k, eps_uk, eps_ud, eps0, beta)
    return float(Mc + Ms)


# ════════════════════════════════════════════════════════════════════════════
# 9. SOLVEUR ELU — (N, M) → (ε₀, β)
# ════════════════════════════════════════════════════════════════════════════

@func
def calculer_N_M_ELU(b,h,bs,hs,gs,bi,hi,gi,asup,ainf,esup,einf, fck,fcd,fyd,k,eps_uk,eps_ud,eps0,beta):
    # Calcul pour la partie béton (Ic et Iv pour les trois valeurs N, My, Mz)
    points = section_I(b,h,bs,hs,gs,bi,hi,gi)
    
    N = N_I_ELU_pararect(b, h, bs, hs, gs, bi, hi, gi,
                      asup, ainf, esup, einf,
                      fck, fcd, fyd, k, eps_uk, eps_ud, eps0, beta)
    M = M_I_ELU_pararect(b, h, bs, hs, gs, bi, hi, gi,
                      asup, ainf, esup, einf,
                      fck, fcd, fyd, k, eps_uk, eps_ud, eps0, beta)
    
    
    return N, M




# ════════════════════════════════════════════════════════════════════════════
# HELPERS INTERNES
# ════════════════════════════════════════════════════════════════════════════

def _residuals_ELU(ep, b, h, bs, hs, gs, bi, hi, gi,
                   asup, ainf, esup, einf,
                   fck, fcd, fyd, k, eps_uk, eps_ud,
                   Nobj, Mobj):
    """Résidu [N−Nobj, M−Mobj] pour un vecteur ep = [ε₀, β]."""
    N, M = calculer_N_M_ELU(
        b, h, bs, hs, gs, bi, hi, gi,
        asup, ainf, esup, einf,
        fck, fcd, fyd, k, eps_uk, eps_ud,
        ep[0], ep[1]
    )
    return np.array([N - Nobj, M - Mobj])


def _jacobian_centre(ep, resid_fn, h_fd=1e-6):
    """
    Jacobien 2×2 par différences finies CENTRÉES O(h²).
    Plus précis que l'unilatéral → meilleure convergence près des kinks.
    """
    J = np.empty((2, 2))
    for i in range(2):
        dh      = np.zeros(2); dh[i] = h_fd
        J[:, i] = (resid_fn(ep + dh) - resid_fn(ep - dh)) / (2.0 * h_fd)
    return J


def _grille_x0(Nobj, Mobj, h):
    """
    Génère 8 points de départ couvrant les domaines physiques EC2 :
      ε₀ ∈ [−3.5, +3.5] ‰  (traction pure → compression pure)
      β  ∈ [−7/h,  +7/h]   (gradient de déformation sur la hauteur)
    """
    eps0_vals = [-2.0, 0.0, 1.0, 2.5]
    beta_vals = [-3.5/h, 0.0, 3.5/h]
    # Point centré sur la sollicitation
    beta_est  = 3.5 / h if Mobj > 0 else -3.5 / h
    eps0_est  = 0.5
    points = [(eps0_est, beta_est)]
    for e in eps0_vals:
        for b_ in beta_vals:
            points.append((e, b_))
    return [np.array([e, b_]) for e, b_ in points]


# ════════════════════════════════════════════════════════════════════════════
# MÉTHODE 1 : root / hybr  (Powell — rapide si régulier)
# ════════════════════════════════════════════════════════════════════════════

def _solve_hybr(resid_fn, x0):
    """root/hybr avec jacobien centré."""
    try:
        sol = root(
            resid_fn, x0,
            jac=lambda ep: _jacobian_centre(ep, resid_fn),
            method="hybr",
            tol=1e-6,
            options={"maxfev": 600}
        )
        if sol.success and np.max(np.abs(sol.fun)) < 1e-3:
            return sol.x
    except Exception:
        pass
    return None


# ════════════════════════════════════════════════════════════════════════════
# MÉTHODE 2 : least_squares / trf  (trust-region — robuste aux kinks)
# ════════════════════════════════════════════════════════════════════════════

def _solve_trf(resid_fn, x0):
    """least_squares/trf — tolère les discontinuités de dérivée."""
    try:
        sol = least_squares(
            resid_fn, x0,
            jac=lambda ep: _jacobian_centre(ep, resid_fn),
            method='trf',
            xtol=1e-6, ftol=1e-6, gtol=1e-6,
            max_nfev=400
        )
        #if np.max(np.abs(sol.fun)) < 1e-3:
        #    return sol.x
    except Exception:
        pass
    return None


# ════════════════════════════════════════════════════════════════════════════
# MÉTHODE 3 : Multi-start  (8 points de départ → hybr puis trf)
# ════════════════════════════════════════════════════════════════════════════

def _solve_multistart(resid_fn, Nobj, Mobj, h):
    """
    Essaie plusieurs points de départ couvrant tout le domaine EC2.
    Pour chaque point : hybr d'abord, puis trf si échec.
    Retourne la meilleure solution (résidu minimal).
    """
    best_x   = None
    best_err = np.inf

    for x0 in _grille_x0(Nobj, Mobj, h):
        # essai hybr
        x = _solve_hybr(resid_fn, x0)
        if x is None:
            x = _solve_trf(resid_fn, x0)
        if x is not None:
            err = float(np.max(np.abs(resid_fn(x))))
            if err < best_err:
                best_err = err
                best_x   = x

    if best_err < 1e-2:
        return best_x
    return None


# ════════════════════════════════════════════════════════════════════════════
# MÉTHODE 4 : Bissection 1D  (filet de sécurité ultime)
# ════════════════════════════════════════════════════════════════════════════

def _solve_bissection(resid_fn, Nobj, Mobj, h):
    """
    Approche séquentielle si toutes les méthodes 2D ont échoué.

    Étape A — trouver β tel que l'équilibre en M est satisfait à ε₀ fixé.
    Étape B — ajuster ε₀ pour l'équilibre en N.

    Hypothèse : β = (ε_cu2 − ε_inf) / h  est monotone en ε₀.
    Fonctionne sur la quasi-totalité des cas EC2 physiquement admissibles.
    """
    from scipy.optimize import brentq as _brentq

    eps_cu2_val = 3.5    # valeur typique C20/C50

    # ── Étape A : chercher β pour un ε₀ donné ────────────────────────────
    def beta_from_eps0(eps0):
        """β tel que le moment calculé = Mobj (à ε₀ fixé)."""
        def fM(beta):
            r = resid_fn(np.array([eps0, beta]))
            return float(r[1])   # M − Mobj
        try:
            return _brentq(fM, -eps_cu2_val / h, eps_cu2_val / h,
                           xtol=1e-8, maxiter=100)
        except ValueError:
            return None

    # ── Étape B : chercher ε₀ tel que l'effort N est satisfait ───────────
    def fN(eps0):
        beta = beta_from_eps0(eps0)
        if beta is None:
            return np.nan
        r = resid_fn(np.array([eps0, beta]))
        return float(r[0])   # N − Nobj

    # Balayage pour trouver un encadrement de ε₀
    eps0_grid = np.linspace(-3.5, 3.5, 30)
    fN_vals   = [fN(e) for e in eps0_grid]

    # Chercher un changement de signe
    for i in range(len(fN_vals) - 1):
        f1, f2 = fN_vals[i], fN_vals[i+1]
        if np.isnan(f1) or np.isnan(f2):
            continue
        if f1 * f2 < 0:
            try:
                eps0_sol = _brentq(fN, eps0_grid[i], eps0_grid[i+1],
                                   xtol=1e-7, maxiter=80)
                beta_sol = beta_from_eps0(eps0_sol)
                if beta_sol is not None:
                    x = np.array([eps0_sol, beta_sol])
                    if np.max(np.abs(resid_fn(x))) < 1e-2:
                        return x
            except Exception:
                continue

    return None


# ════════════════════════════════════════════════════════════════════════════
# SOLVEUR PRINCIPAL — Cascade des 4 méthodes
# ════════════════════════════════════════════════════════════════════════════

@func
def solve_I_ELU_pararect(
    b, h, bs, hs, gs, bi, hi, gi,
    asup, ainf, esup, einf,
    fck, fcd, fyd, k, eps_uk, eps_ud,
    Nobj, Mobj
):
    """
    Résout (N = Nobj, M = Mobj) → (ε₀, β) — ELU parabole-rectangle.

    Cascade de robustesse :
      1. root/hybr          (le plus rapide)
      2. least_squares/trf  (robuste aux kinks β/ε₀)
      3. Multi-start        (8 points de départ différents)
      4. Bissection 1D      (filet de sécurité ultime)

    Retourne [ε₀, β] ou [nan, nan] si aucune méthode ne converge.
    """

    # Fermeture sur les paramètres de la section
    def resid(ep):
        return _residuals_ELU(
            ep, b, h, bs, hs, gs, bi, hi, gi,
            asup, ainf, esup, einf,
            fck, fcd, fyd, k, eps_uk, eps_ud,
            Nobj, Mobj
        )

    x0 = np.array([0.01, 0.1])   # point de départ par défaut

    # ── Méthode 1 : hybr ─────────────────────────────────────────────────
    x = _solve_hybr(resid, x0)
    if x is not None:
        return x

    # ── Méthode 2 : trf ──────────────────────────────────────────────────
    x = _solve_trf(resid, x0)
    if x is not None:
        return x

    # ── Méthode 3 : multi-start ───────────────────────────────────────────
    x = _solve_multistart(resid, Nobj, Mobj, h)
    if x is not None:
        return x

    # ── Méthode 4 : bissection 1D ─────────────────────────────────────────
    x = _solve_bissection(resid, Nobj, Mobj, h)
    if x is not None:
        return x

    # Aucune convergence
    return np.array([np.nan, np.nan])


# ════════════════════════════════════════════════════════════════════════════
# 10. Résultats ELU — Calculs contrainte
# ════════════════════════════════════════════════════════════════════════════   
@func
def resultats_I_ELU_pararect(b,h,bs,hs,gs,bi,hi,gi,asup,ainf,esup,einf, fck,fcd,fyd,k,eps_uk,eps_ud,eps0,beta):
    # Résultats de calcul sur une section soumise à un champ de déformation eps = eps0+beta.y
    epscsup = eps0+(h-Nc_Gy_ELS(b, h, bs, hs,gs,bi,hi,gi))*beta
    epscinf = eps0-Nc_Gy_ELS(b, h, bs, hs,gs,bi,hi,gi)*beta
    epsssup = eps0+(h-Nc_Gy_ELS(b, h, bs, hs,gs,bi,hi,gi)-esup)*beta
    epssinf = eps0-(Nc_Gy_ELS(b, h, bs, hs,gs,bi,hi,gi)-einf)*beta
    sigmacsup = float(sigma_c_pararect1(float(fck), float(fcd), epscsup))
    sigmacinf = float(sigma_c_pararect1(float(fck), float(fcd), epscinf))

    sigmassup = float(sigma_s_palier(fyd, k, eps_uk, eps_ud, epsssup))
    sigmasinf = float(sigma_s_palier(fyd, k, eps_uk, eps_ud, epssinf))

    N, M =calculer_N_M_ELU(b,h,bs,hs,gs,bi,hi,gi,asup,ainf,esup,einf, fck,fcd,fyd,k,eps_uk,eps_ud,eps0,beta)

    ycomp = lambda y: np.where((eps0 + beta * y) > 0, 1.0, 0.0)
    hcomp = scipy.integrate.quad(ycomp, -Nc_Gy_ELS(b, h, bs, hs,gs,bi,hi,gi), h-Nc_Gy_ELS(b, h, bs, hs,gs,bi,hi,gi))[0]

    if hcomp==h:
        etat="entierement comprimé"
    elif hcomp>0:
        etat="partiellement tendu"
    else:
        etat="entierement tendu"

    res = {}
    res["EPS0"]=eps0
    res["BETA"]=beta
    res["H_compr"]=hcomp
    res["Etat"]= etat
    res["EPS_C_SUP"]=epscsup
    res["EPS_C_INF"]=epscinf
    res["EPS_S_SUP"]=epsssup
    res["EPS_S_INF"]=epssinf
    res["SIG_C_SUP"]=sigmacsup
    res["SIG_C_INF"]=sigmacinf
    res["SIG_S_SUP"]=sigmassup
    res["SIG_S_INF"]=sigmasinf
    res["N"]=N
    res["M"]=M

    return res

@func
def e_resultats_I_ELU_pararect(b,h,bs,hs,gs,bi,hi,gi,asup,ainf,esup,einf, fck,fcd,fyd,k,eps_uk,eps_ud,eps0,beta,resultats):
    resultats_Excel = []
    resultats_tout = resultats_I_ELU_pararect(b,h,bs,hs,gs,bi,hi,gi,asup,ainf,esup,einf, fck,fcd,fyd,k,eps_uk,eps_ud,eps0,beta)                     
    resultats_list = resultats.split(',')
    for r in resultats_list:       
        resultats_Excel.append(resultats_tout[r])
    return resultats_Excel   # résultats en ligne
# ════════════════════════════════════════════════════════════════════════════
# 11. MOMENTS DE RÉFÉRENCE ELS
# ════════════════════════════════════════════════════════════════════════════

def _eps_beta_AB_ELS(yG, h, einf, n, sb, syt):
    """Plan de déformation pivotant sur fibre tendue + acier inf à syt."""
    beta = (n * sb / 200.0 - syt / 200.0) / (h - einf)
    eps0 = syt / 200.0 + (yG - einf) * beta
    return eps0, beta


@func
def ELS_I_MserA(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Moment de service par rapport à l'acier inf."""
    yG = _yG(b, h, bs, hs, gs, bi, hi, gi)
    return Ms + Ns * (yG - einf)


@func
def ELS_I_MserB(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Moment de service par rapport à l'acier sup."""
    yG = _yG(b, h, bs, hs, gs, bi, hi, gi)
    return -Ms + Ns * (h - yG - esup)


@func
def ELS_I_MAB(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Moment de référence MAB (pivot acier inf = syt)."""
    yG       = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts      = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    eps0, bt = _eps_beta_AB_ELS(yG, h, einf, n, sb, syt)
    Nc, Mc   = _NM_beton_ELS(pts, n, eps0, bt)
    return float(Mc + Nc * (yG - einf))

@func
def ELS_I_NAB(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Moment de référence MAB (pivot acier inf = syt)."""
    yG       = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts      = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    eps0, bt = _eps_beta_AB_ELS(yG, h, einf, n, sb, syt)
    Nc, _  = _NM_beton_ELS(pts, n, eps0, bt)
    return float(Nc) 



@func
def ELS_I_MAB_p(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Moment de référence MAB' (par rapport à acier sup)."""
    yG       = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts      = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    eps0, bt = _eps_beta_AB_ELS(yG, h, einf, n, sb, syt)
    Nc, Mc   = _NM_beton_ELS(pts, n, eps0, bt)
    return float(Nc * (-Mc / Nc + h - yG - esup))


@func
def ELS_I_MBO(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Moment de référence MBO (pivot fibre inf = 0)."""
    yG     = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts    = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    beta   = n * sb / 200.0 / h
    eps0   = yG * beta
    Nc, Mc = _NM_beton_ELS(pts, n, eps0, beta)
    return float(Nc * (-Mc / Nc + h - yG - esup))


@func
def ELS_I_MBMAX(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Moment de référence MBMAX (section entièrement comprimée)."""
    yG     = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts    = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    Nc, Mc = _NM_beton_ELS(pts, n, n * sb / 200.0, 0.0)
    return float(Nc * (-Mc / Nc + h - yG - esup))


# ════════════════════════════════════════════════════════════════════════════
# 12. CALCUL DES ARMATURES ELS — MÉTHODE DES CAS
# ════════════════════════════════════════════════════════════════════════════

def _brentq(f, a, b, tol=1e-12):
    """Encapsule root_scalar Brent — retourne sol.root ou valeur de repli b."""
    try:
        return root_scalar(f, bracket=[a, b], method='brentq', xtol=tol).root
    except ValueError:
        return b


@func
def ELS_I_As_t(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Cas 0 — section entièrement tendue : Asup."""
    MserA = ELS_I_MserA(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    return abs(MserA) / (h - esup - einf) / (-syt) * 1e4


@func
def ELS_I_Ai_t(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Cas 0 — section entièrement tendue : Ainf."""
    MserA = ELS_I_MserA(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    return Ns / syt * 1e4 - abs(MserA) / (h - esup - einf) / (-syt) * 1e4


# ── Cas 1 : Asup = 0 ─────────────────────────────────────────────────────

@func
def solve_I_ELS_c1(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Résout y (profondeur relative) pour le cas 1."""
    MserA = ELS_I_MserA(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    pts   = _pts_I(b, h, bs, hs, gs, bi, hi, gi)

    def f(y):
        y     = np.clip(y, 1e-6, 1.0 - 1e-6)
        sig_b = syt / n * (1.0 - 1.0 / y)
        eps0_, bt_ = _eps_beta_AB_ELS(
            _yG(b, h, bs, hs, gs, bi, hi, gi), h, einf, n, sig_b, syt)
        Nc_, Mc_ = _NM_beton_ELS(pts, n, eps0_, bt_)
        yG_ = _yG(b, h, bs, hs, gs, bi, hi, gi)
        return MserA - (Mc_ + Nc_ * (yG_ - einf))

    return _brentq(f, 1e-5, 1.0 - 1e-5)


@func
def ELS_I_Ai_1(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Cas 1 — Asup = 0 : Ainf."""
    y        = solve_I_ELS_c1(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    yG       = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts      = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    sig_b    = syt / n * (1.0 - 1.0 / y)
    eps0, bt = _eps_beta_AB_ELS(yG, h, einf, n, sig_b, syt)
    Nc, _    = _NM_beton_ELS(pts, n, eps0, bt)
    return float((Ns - Nc) / syt * 1e4)


# ── Cas 2 : Asup > 0, Ainf > 0 ───────────────────────────────────────────

@func
def ELS_I_sAs_2(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Contrainte Asup — cas 2."""
    alpha = n*sb/(n*sb - syt)
    sig_s = n*sb*(1-esup/(h-einf)/alpha)
    return min(syc, sig_s)


@func
def ELS_I_As_2(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Cas 2 — Asup."""
    MserA = ELS_I_MserA(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    MAB   = ELS_I_MAB(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    sAs   = ELS_I_sAs_2(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    return float(1e4 * (MserA - MAB) / (h - esup - einf) / sAs)


@func
def ELS_I_Ai_2(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Cas 2 — Ainf."""
    #MAB    = ELS_I_MAB(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    NAB    = ELS_I_NAB(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    sAs    = ELS_I_sAs_2(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    As_2   = ELS_I_As_2(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    
     
    return float(1e4 / (-syt) * (NAB+ As_2*1e-4 * sAs - Ns))


# ── Cas 22 : Asup > 0, Ainf = 0 ──────────────────────────────────────────

@func
def solve_I_ELS_c22(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Résout y (profondeur axe neutre) pour le cas 22."""
    MserB = ELS_I_MserB(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    yG    = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts   = _pts_I(b, h, bs, hs, gs, bi, hi, gi)

    def f(y):
        y_s  = max(y, 1e-4)
        s_var = sb * n * (y_s - h + einf) / y_s
        eps0_, bt_ = _eps_beta_AB_ELS(yG, h, einf, n, sb, s_var)
        Nc_, Mc_   = _NM_beton_ELS(pts, n, eps0_, bt_)
        return MserB - (-Mc_ + Nc_ * (h - yG - esup))

    return _brentq(f, 1e-4, h * 5.0)


@func
def ELS_I_sAs_22(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    y = solve_I_ELS_c22(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    return min(syc, n * sb * (1.0 - esup / y))


@func
def ELS_I_As_22(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Cas 22 — Asup, Ainf = 0."""
    y        = solve_I_ELS_c22(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    yG       = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts      = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    s_var    = sb * n * (1.0 - (h - einf) / y)
    eps0, bt = _eps_beta_AB_ELS(yG, h, einf, n, sb, s_var)
    Nc, _    = _NM_beton_ELS(pts, n, eps0, bt)
    sAs      = ELS_I_sAs_22(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    return float((Ns - Nc) / sAs * 1e4)


# ── Cas 3 : pivot fibre inf (Ainf = 0) ───────────────────────────────────

@func
def ELS_I_MBO3(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    beta= (n * sb / 200 - syt / 200) / (h - einf)
    epo= syt / 200 + (Nc_Gy_ELS(b, h, bs, hs, gs, bi, hi, gi) - einf)*beta     
    Mc_value = Mc_I_ELS(b, h, bs, hs, gs, bi, hi, gi, n, epo, beta)
    Nc_value = Nc_I_ELS(b, h, bs, hs, gs, bi, hi, gi, n, epo, beta)
    Nc_Gy_value = Nc_Gy_ELS(b, h, bs, hs, gs, bi, hi, gi)

    return -Mc_value + Nc_value * (h - Nc_Gy_value - esup)
    
@func    
def solve_I_ELS_c3(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
  
    y0 = 4*h  # valeur initiale 
    def fN(y):
        return ELS_I_MserB(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms) - ELS_I_MBO3(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, sb * n * (y - h + einf) / y, syc, Ns, Ms)
    y_solution = fsolve(fN, y0)
    
    return y_solution[0]  


@func
def ELS_I_sAs_3(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    y = solve_I_ELS_c3(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    return min(syc, n * sb * (1.0 - esup / y))


@func
def ELS_I_As_3(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Cas 3."""
    y        = solve_I_ELS_c3(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    yG       = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts      = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    s_var    = sb * n * (1.0 - (h - einf) / y)
    eps0, bt = _eps_beta_AB_ELS(yG, h, einf, n, sb, s_var)
    Nc, _    = _NM_beton_ELS(pts, n, eps0, bt)
    sAs      = ELS_I_sAs_3(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    return float((Ns - Nc) / sAs * 1e4)


# ── Cas 4 : section entièrement comprimée ────────────────────────────────

@func
def ELS_I_Ai_4(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    MserB = ELS_I_MserB(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    MBMax = ELS_I_MBMAX(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    return float((MserB - MBMax) / n / sb / (h - esup - einf) * 1e4)


@func
def ELS_I_As_4(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    Ai_4 = ELS_I_Ai_4(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    pts  = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    Nc, _ = _NM_beton_ELS(pts, n, n * sb / 200.0, 0.0)
    return float((Ns - float(Nc)) / n / sb * 1e4 - Ai_4)
# ── Sélection automatique du cas ELS ─────────────────────────────────────

@func
def section_area(b, h, bs, hs, gs, bi, hi, gi):
    """
    Calcule l'aire totale de la section en I avec goussets.
    Calcul analytique exact (O(1)), sans polygone_integrate ni sommets.
    """
    # 1. Aire de la table inférieure (Rectangle)
    A_table_inf = bi * hi
    
    # 2. Aire du gousset inférieur (Trapèze)
    A_gousset_inf = 0.5 * (bi + b) * gi
    
    # 3. Aire de l'âme (Rectangle)
    # Sa hauteur est la hauteur totale moins les tables et les goussets
    h_ame = h - hi - gi - hs - gs
    A_ame = b * h_ame
    
    # 4. Aire du gousset supérieur (Trapèze)
    A_gousset_sup = 0.5 * (b + bs) * gs
    
    # 5. Aire de la table supérieure (Rectangle)
    A_table_sup = bs * hs
    
    # Somme de toutes les sous-surfaces
    aire_totale = A_table_inf + A_gousset_inf + A_ame + A_gousset_sup + A_table_sup
    
    return float(aire_totale)


@func
def ELS_I_As_M(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Sélection automatique du cas ELS et calcul Asup."""
    yG = _yG(b, h, bs, hs, gs, bi, hi, gi)
    if Ns < 0 and Ms / (Ns - 1e-15) >= -(yG - einf):
        return ELS_I_As_t(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)

    MserA = ELS_I_MserA(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    MAB   = ELS_I_MAB(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    if MserA <= MAB:
        return 0.0

    MserB = ELS_I_MserB(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    MABp  = ELS_I_MAB_p(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    if MserB <= MABp:
        return ELS_I_As_2(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)

    MBO = ELS_I_MBO(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    if MserB <= MBO:
        return ELS_I_As_22(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)

    MBMax = ELS_I_MBMAX(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    if MserB <= MBMax:
        return ELS_I_As_3(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)

    return ELS_I_As_4(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)


@func
def ELS_I_Ai_M(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    """Sélection automatique du cas ELS et calcul Ainf."""
    yG = _yG(b, h, bs, hs, gs, bi, hi, gi)
    if Ns < 0 and Ms / (Ns - 1e-15) >= -(yG - einf):
        return ELS_I_Ai_t(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)

    MserA = ELS_I_MserA(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    MAB   = ELS_I_MAB(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    if MserA <= MAB:
        return ELS_I_Ai_1(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)

    MserB = ELS_I_MserB(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    MABp  = ELS_I_MAB_p(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    if MserB <= MABp:
        return ELS_I_Ai_2(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)

    MBO = ELS_I_MBO(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    if MserB <= MBO:
        return 0.0

    MBMax = ELS_I_MBMAX(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    if MserB <= MBMax:
        return 0.0

    return ELS_I_Ai_4(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)


@func
def ELS_I_As_Max(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    if Ms > 0:
        return ELS_I_As_M(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    return ELS_I_Ai_M(b, h, bi, hi, gi, bs, hs, gs, einf, esup, n, sb, syt, syc, Ns, -Ms)


@func
def ELS_I_Ai_Max(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    if Ms > 0:
        return ELS_I_Ai_M(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    return ELS_I_As_M(b, h, bi, hi, gi, bs, hs, gs, einf, esup, n, sb, syt, syc, Ns, -Ms)


def _check_4pct(As, Ai, area):
    """Retourne un message si la section totale dépasse 4 %."""
    if (As + Ai) > 400.0 * area:
        return "Section totale dépasse 4 %"
    return None


@func
def ELS_I_As(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    As = ELS_I_As_Max(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    Ai = ELS_I_Ai_Max(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    msg = _check_4pct(As, Ai, section_area(b, h, bs, hs, gs, bi, hi, gi))
    if msg:   return msg
    return max(0.0, As)


@func
def ELS_I_Ai(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    As = ELS_I_As_Max(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    Ai = ELS_I_Ai_Max(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms)
    msg = _check_4pct(As, Ai, section_area(b, h, bs, hs, gs, bi, hi, gi))
    if msg:   return msg
    return max(0.0, Ai)


@func
def ELS_I_A(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms):
    return (ELS_I_As(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms),
            ELS_I_Ai(b, h, bs, hs, gs, bi, hi, gi, esup, einf, n, sb, syt, syc, Ns, Ms))


# ════════════════════════════════════════════════════════════════════════════
# 13. MOMENTS DE RÉFÉRENCE ELU
# ════════════════════════════════════════════════════════════════════════════

def _beta_epo_ELU(yG, h, einf, eps_bot, eps_top, pts, fck, fcd):
    """Plan de déformation ELU générique → (eps0, beta, Nc, Mc)."""
    beta   = (eps_top + eps_bot) / (h - einf)
    eps0   = -eps_bot + (yG - einf) * beta
    Nc, Mc = _NM_beton_ELU(pts, fck, fcd, eps0, beta)
    return eps0, beta, float(Nc), float(Mc)


@func
def ELU_I_MuA(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    """Moment appliqué / acier sup (bras de levier h−yG−esup)."""
    yG = _yG(b, h, bs, hs, gs, bi, hi, gi)
    return float(-Mu + Nu * (h - yG - esup))


@func
def ELU_I_MuA1(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    """Moment appliqué / acier inf (bras de levier yG−einf)."""
    yG = _yG(b, h, bs, hs, gs, bi, hi, gi)
    return float(Mu + Nu * (yG - einf))


@func
def ELU_I_MAB(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    """Moment de référence MAB (pivot eps_ud à l'acier inf)."""
    yG  = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    _, _, Nc, Mc = _beta_epo_ELU(yG, h, einf, eps_ud, eps_cu2(fck), pts, fck, fcd)
    return float(Mc + Nc * (yG - einf))


@func
def ELU_I_ME(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    """Moment de référence ME (pivot fyd/Es à l'acier inf)."""
    yG  = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    _, _, Nc, Mc = _beta_epo_ELU(yG, h, einf, fyd / 200.0, eps_cu2(fck), pts, fck, fcd)
    return float(Mc + Nc * (yG - einf))


@func
def ELU_I_ME_p(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    """ME' — par rapport à acier sup."""
    yG  = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    _, _, Nc, Mc = _beta_epo_ELU(yG, h, einf, fyd / 200.0, eps_cu2(fck), pts, fck, fcd)
    return float(-Mc + Nc * (h - yG - esup))


@func
def ELU_I_MBC_p(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    """MBC' — pivot axe neutre = fibre inf."""
    yG  = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    beta   = eps_cu2(fck) / h
    eps0   = yG * beta
    Nc, Mc = _NM_beton_ELU(pts, fck, fcd, eps0, beta)
    return float(-Mc + Nc * (h - yG - esup))


@func
def ELU_I_M2_p(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    """M2% — section entièrement comprimée (eps = eps_c2, beta = 0)."""
    yG  = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    Nc, Mc = _NM_beton_ELU(pts, fck, fcd, float(eps_c2(fck)), 0.0)
    return float(-Mc + Nc * (h - yG - esup))


# ════════════════════════════════════════════════════════════════════════════
# 14. CALCUL DES ARMATURES ELU — MÉTHODE DES DOMAINES
# ════════════════════════════════════════════════════════════════════════════

def _solve_ELU_brentq(b, h, bs, hs, gs, bi, hi, gi, fck, fcd,
                       target_M, bras, eps_c_ref, lo, hi_b):
    """
    Résout l'équilibre en moment (target_M = Mc + Nc × bras)
    en cherchant la profondeur d'axe neutre x dans [lo, hi_b].
    Retourne x, eps0, beta, Nc, Mc.
    """
    yG  = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts = _pts_I(b, h, bs, hs, gs, bi, hi, gi)

    def f(x):
        x    = max(x, 1e-6)
        beta = eps_c_ref / x
        eps0 = eps_c_ref - (h - yG) * beta
        Nc_, Mc_ = _NM_beton_ELU(pts, fck, fcd, eps0, beta)
        return target_M - Mc_ - Nc_ * bras

    x = _brentq(f, lo, hi_b)
    beta = eps_c_ref / x
    eps0 = eps_c_ref - (h - yG) * beta
    Nc, Mc = _NM_beton_ELU(pts, fck, fcd, eps0, beta)
    return x, eps0, beta, float(Nc), float(Mc)


# ── Cas tendu (tout acier) ────────────────────────────────────────────────

@func
def ELU_I_As_t(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    sig_s=sigma_s_palier(fyd, k, eps_uk,eps_ud, eps_ud) # allongement max
    MuA1 = abs(ELU_I_MuA1(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu))
    return float(MuA1 / (h - esup - einf) / sig_s * 1e4)


@func
def ELU_I_Ai_t(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    sig_s=sigma_s_palier(fyd, k, eps_uk,eps_ud, eps_ud) # allongement max
    MuA1 = abs(ELU_I_MuA1(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu))
    return float(Nu / (-sig_s) * 1e4 - MuA1 / (h - esup - einf) / sig_s * 1e4)


# ── Domaine 1 : eps_inf = eps_ud (très grande traction) ──────────────────

@func
def solve_I_ELU_c1(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    MuA1 = ELU_I_MuA1(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    yG   = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts  = _pts_I(b, h, bs, hs, gs, bi, hi, gi)

    def f(alp):
        alp  = np.clip(alp, 1e-6, 0.999 * (h - einf))
        denom = h - einf - alp
        beta = eps_ud / denom
        eps0 = -eps_ud + (yG - einf) * beta
        Nc_, Mc_ = _NM_beton_ELU(pts, fck, fcd, eps0, beta)
        return MuA1 - Mc_ - Nc_ * (yG - einf)

    return _brentq(f, 1e-6, 0.999 * (h - einf))


@func
def ELU_I_Ai_1(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    alp  = solve_I_ELU_c1(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    yG   = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts  = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    beta = eps_ud / (h - einf - alp)
    eps0 = -eps_ud + (yG - einf) * beta
    Nc, _   = _NM_beton_ELU(pts, fck, fcd, eps0, beta)
    sig_s   = float(sigma_s_palier(fyd, k, eps_uk, eps_ud, eps_ud))
    return float(-(Nu - Nc) / sig_s * 1e4)


# ── Domaine 2 : pivot eps_cu2 en fibre sup ───────────────────────────────

@func
def solve_I_ELU_c2(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    MuA1 = ELU_I_MuA1(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    yG   = _yG(b, h, bs, hs, gs, bi, hi, gi)
    bras = yG - einf
    return _solve_ELU_brentq(b, h, bs, hs, gs, bi, hi, gi, fck, fcd,
                              MuA1, bras, float(eps_cu2(fck)), 1e-4, 5.0*h)[0]


@func
def ELU_I_Ai_2(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    alp2 = solve_I_ELU_c2(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    yG   = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts  = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    ecu  = float(eps_cu2(fck))
    beta = ecu / alp2
    eps0 = ecu - (h - yG) * beta
    Nc, _ = _NM_beton_ELU(pts, fck, fcd, eps0, beta)
    sig_s = float(sigma_s_palier(fyd, k, eps_uk, eps_ud, eps0 + beta * (-yG + einf)))
    return float((Nu - Nc) / sig_s * 1e4)


# ── Domaine 2.1 : Asup > 0, Ainf > 0 (pivot ME) ─────────────────────────

@func
def ELU_I_As_21(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    MuA1  = ELU_I_MuA1(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    ME    = ELU_I_ME(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    ecu   = float(eps_cu2(fck))
    eps_sd = fyd / 200000 * 1000.0
    bras_x = einf * (1.0 + 200.0 * ecu / fyd) * fyd / 200.0 / ecu / (h - einf)
    sig_s  = float(sigma_s_palier(fyd, k, eps_uk, eps_ud, ecu * (1.0 - bras_x)))
    return float(1e4 * (MuA1 - ME) / (h - esup - einf) / sig_s)


@func
def ELU_I_Ai_21(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    yG    = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts   = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    MuA1  = ELU_I_MuA1(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    ME    = ELU_I_ME(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    ecu   = float(eps_cu2(fck))
    beta  = (ecu + fyd / 200.0) / (h - einf)
    eps0  = -fyd / 200.0 + (yG - einf) * beta
    Nc, _ = _NM_beton_ELU(pts, fck, fcd, eps0, beta)
    return float(1e4 / fyd * (-Nu + Nc + (MuA1 - ME) / (h - esup - einf)))


# ── Domaine 2.2 : Asup > 0, Ainf = 0 ────────────────────────────────────

@func
def solve_I_ELU_c22(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    ELU_MA = ELU_I_MuA(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    Nc_gy = Nc_Gy_ELS(b, h, bs, hs, gs, bi, hi, gi)
 
    
    x0 = 0.8*h  # valeur initiale 
    def fN(x):
        Mc_I_ELU= Mc_I_ELU_pararect(b, h, bs, hs, gs, bi, hi, gi, fck, fcd, eps_cu2(fck) - (h - Nc_gy) * eps_cu2(fck) / x, eps_cu2(fck) / x)
        Nc_I_ELU = Nc_I_ELU_pararect(b, h, bs, hs, gs, bi, hi, gi, fck, fcd, eps_cu2(fck) - (h - Nc_gy) * eps_cu2(fck) / x, eps_cu2(fck) / x)
        return ELU_MA + Mc_I_ELU - Nc_I_ELU * (h - Nc_gy - esup)

    x_solution = fsolve(fN, x0)
    return x_solution[0]


@func
def ELU_I_As_22(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    alp  = solve_I_ELU_c22(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    yG   = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts  = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    ecu  = float(eps_cu2(fck))
    beta = ecu / alp
    eps0 = ecu - (h - yG) * beta
    Nc, _ = _NM_beton_ELU(pts, fck, fcd, eps0, beta)
    sig_s = float(sigma_s_palier(fyd, k, eps_uk, eps_ud, eps0 + beta * (h - yG - esup)))
    return float((Nu - Nc) / sig_s * 1e4)


# ── Domaine 3 : M2% > MuA' (section très comprimée) ─────────────────────


@func
def solve_I_ELU_c3(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    ELU_MA = ELU_I_MuA(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    Nc_gy = Nc_Gy_ELS(b, h, bs, hs, gs, bi, hi, gi)
 
    
    x0 = 0.8*h  # valeur initiale 
    def fN(x):
        Mc_I_ELU= Mc_I_ELU_pararect(b, h, bs, hs, gs, bi, hi, gi, fck, fcd, eps_cu2(fck) - (h - Nc_gy) * eps_cu2(fck) / x, eps_cu2(fck) / x)
        Nc_I_ELU = Nc_I_ELU_pararect(b, h, bs, hs, gs, bi, hi, gi, fck, fcd, eps_cu2(fck) - (h - Nc_gy) * eps_cu2(fck) / x, eps_cu2(fck) / x)
        return ELU_MA + Mc_I_ELU - Nc_I_ELU * (h - Nc_gy - esup)

    x_solution = fsolve(fN, x0)
    return x_solution[0] 

@func
def ELU_I_As_3(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    alp  = solve_I_ELU_c3(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    yG   = _yG(b, h, bs, hs, gs, bi, hi, gi)
    pts  = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    ecu  = float(eps_cu2(fck))
    beta = ecu / alp
    eps0 = ecu - (h - yG) * beta
    Nc, _ = _NM_beton_ELU(pts, fck, fcd, eps0, beta)
    sig_s = float(sigma_s_palier(fyd, k, eps_uk, eps_ud, ecu * (alp - esup) / alp))
    return float((Nu - Nc) / sig_s * 1e4)

@func
def ELU_I_Ai_4(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    MuA   = ELU_I_MuA(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    M2p   = ELU_I_M2_p(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    sig_s = float(sigma_s_palier(fyd, k, eps_uk, eps_ud, float(eps_c2(fck))))
    return float(1e4 * (MuA - M2p) / (h - esup - einf) / sig_s)


@func
def ELU_I_As_4(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    pts   = _pts_I(b, h, bs, hs, gs, bi, hi, gi)
    Nc, _ = _NM_beton_ELU(pts, fck, fcd, float(eps_c2(fck)), 0.0)
    MuA   = ELU_I_MuA(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    M2p   = ELU_I_M2_p(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    sig_s = float(sigma_s_palier(fyd, k, eps_uk, eps_ud, float(eps_c2(fck))))
    return float(1e4 / sig_s * (Nu - Nc - (MuA - M2p) / (h - esup - einf)))

# ── Sélection automatique du domaine ELU ─────────────────────────────────

@func
def ELU_I_As_M(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    """Asup — sélection automatique du domaine."""
    yG   = _yG(b, h, bs, hs, gs, bi, hi, gi)
    args = (b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)

    if Nu < 0 and Mu / (Nu - 1e-15) >= -(yG - einf):
        return ELU_I_As_t(*args)

    MuA1 = ELU_I_MuA1(*args);  ME   = ELU_I_ME(*args)
    if MuA1 <= ME:
        return 0.0

    MuA  = ELU_I_MuA(*args);   MEp  = ELU_I_ME_p(*args)
    if MuA <= MEp:
        return ELU_I_As_21(*args)

    MBCp = ELU_I_MBC_p(*args)
    if MuA <= MBCp:
        return ELU_I_As_22(*args)

    M2p = ELU_I_M2_p(*args)
    if MuA <= M2p:
        return ELU_I_As_3(*args)

    return ELU_I_As_4(*args)


@func
def ELU_I_Ai_M(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    """Ainf — sélection automatique du domaine."""
    yG   = _yG(b, h, bs, hs, gs, bi, hi, gi)
    args = (b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)

    if Nu < 0 and Mu / (Nu - 1e-15) >= -(yG - einf):
        return ELU_I_Ai_t(*args)

    MuA1 = ELU_I_MuA1(*args);  MAB = ELU_I_MAB(*args)
    if MuA1 <= MAB:
        return ELU_I_Ai_1(*args)

    ME = ELU_I_ME(*args)
    if MuA1 <= ME:
        return ELU_I_Ai_2(*args)

    MuA = ELU_I_MuA(*args);  MEp = ELU_I_ME_p(*args)
    if MuA <= MEp:
        return ELU_I_Ai_21(*args)

    MBCp = ELU_I_MBC_p(*args)
    if MuA <= MBCp:
        return 0.0

    M2p = ELU_I_M2_p(*args)
    if MuA <= M2p:
        return 0.0

    return ELU_I_Ai_4(*args)


@func
def ELU_I_As_Max(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    if Mu > 0:
        return ELU_I_As_M(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    return ELU_I_Ai_M(b, h, bi, hi, gi, bs, hs, gs, einf, esup, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, -Mu)


@func
def ELU_I_Ai_Max(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    if Mu > 0:
        return ELU_I_Ai_M(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    return ELU_I_As_M(b, h, bi, hi, gi, bs, hs, gs, einf, esup, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, -Mu)


@func
def ELU_I_As(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    As = ELU_I_As_Max(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    Ai = ELU_I_Ai_Max(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    msg = _check_4pct(As, Ai, section_area(b, h, bs, hs, gs, bi, hi, gi))
    if msg:   return msg
    return max(0.0, As)


@func
def ELU_I_Ai(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    As = ELU_I_As_Max(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    Ai = ELU_I_Ai_Max(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu)
    msg = _check_4pct(As, Ai, section_area(b, h, bs, hs, gs, bi, hi, gi))
    if msg:   return msg
    return max(0.0, Ai)


@func
def ELU_I_A(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu):
    return (ELU_I_As(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu),
            ELU_I_Ai(b, h, bs, hs, gs, bi, hi, gi, esup, einf, fck, fcd, fyd, k, eps_uk, eps_ud, Nu, Mu))

# Dessin la section
import matplotlib.pyplot as plt
import io
@func
def dessiner_section_I(b, h, bs, hs, gs, bi, hi, gi):
    # 1. Tes calculs géométriques (appels à tes fonctions existantes)
    yg_val = Nc_Gy_ELS(b, h, bs, hs, gs, bi, hi, gi)
    poly_centre = section_I(b, h, bs, hs, gs, bi, hi, gi)

    # 2. Création du graphique avec le backend non-interactif
    fig, ax = plt.subplots(figsize=(5, 7))
    
    x_coords = poly_centre[:, 0]
    y_coords = poly_centre[:, 1]
    
    ax.plot(x_coords, y_coords, color='navy', linewidth=2)
    ax.fill(x_coords, y_coords, color='lightsteelblue', alpha=0.5)
    ax.axhline(0, color='red', linestyle='--', linewidth=1)
    ax.set_aspect('equal')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.set_title(f"Section Recentrée (yg={yg_val:.3f})")

    # 3. CRUCIAL POUR LITE : Convertir la figure en image pour le retour
    # Si tu utilises le système de fonctions personnalisées (UDF) de Lite :
    return fig
@func
def create_triangles(vertices, max_depth=0):
    vertices = np.array(vertices)
    tri = Delaunay(vertices, qhull_options="QJ")
    polygon_path = Path(vertices)
    
    def is_triangle_inside(triangle):
        centroid = np.mean(triangle, axis=0)
        return polygon_path.contains_point(centroid)
    
    def subdivide_triangle(triangle, depth):
        if depth >= max_depth:
            return [triangle]
        
        midpoints = [(triangle[i] + triangle[(i + 1) % 3]) / 2 for i in range(3)]
        new_triangles = [
            [triangle[0], midpoints[0], midpoints[2]],
            [triangle[1], midpoints[0], midpoints[1]],
            [triangle[2], midpoints[1], midpoints[2]],
            [midpoints[0], midpoints[1], midpoints[2]]
        ]
        
        subdivided_triangles = []
        for t in new_triangles:
            subdivided_triangles.extend(subdivide_triangle(np.array(t), depth + 1))
        
        return subdivided_triangles
    
    triangles = []
    for simplex in tri.simplices:
        triangle = vertices[simplex]
        if is_triangle_inside(triangle):
            triangles.extend(subdivide_triangle(triangle, 0))
    
    return triangles
@func
@func
def plot_triangles(b, h, bs, hs, gs, bi, hi, gi, Asup, Ainf, esup, einf, ratio):
    # 1. Calculs géométriques
    points = section_I(b, h, bs, hs, gs, bi, hi, gi)
    tri_pol = create_triangles(points)    
    yg = Nc_Gy_ELS(b, h, bs, hs, gs, bi, hi, gi)
    
    # 2. Création de la figure (Nettoyage des doubles fenêtres)
    fig, ax = plt.subplots()
    #fig, ax = plt.subplots(figsize=(6, 6))

    # 3. Tracé du contour bleu
    x, y = zip(*points)
    x_plot = list(x) + [x[0]]
    y_plot = list(y) + [y[0]]
    ax.plot(x_plot, y_plot, linestyle='-', color='blue', linewidth=1.0)

    # 4. Remplissage par triangles (Gris transparent)
    for triangle in tri_pol:
        t_arr = np.array(triangle)
        ax.fill(t_arr[:, 0], t_arr[:, 1], 
                edgecolor=(0.5, 0.5, 0.5, 0.2), 
                facecolor=(0.6, 0.6, 0.6, 0.5), 
                linewidth=0.1)

    # 5. Dessin des Aciers (Points seuls, sans ligne)
    # On définit un décalage horizontal dynamique pour le texte (5% de la largeur)
    offset = b * 0.1 

    # Acier Supérieur (Vert)
    yas = h - yg - esup
    ax.plot(0, yas, marker='o', color='green', markersize=5, linestyle='None')
    ax.text(offset, yas,  f"Asup", 
            color='black', fontsize=8, va='center')

    # Acier Inférieur (Orange)
    yai = einf - yg 
    ax.plot(0, yai, marker='o', color='orange', markersize=5, linestyle='None')
    ax.text(offset, yai, f"Ainf", 
            color='black', fontsize=8, va='center')

    # 6. Configuration finale de l'affichage
    # set_aspect(1/ratio) pour respecter les unités réelles demandées
    ax.set_aspect(ratio, adjustable='datalim')
    
    ax.set_title('Section prise en compte', fontweight='bold', pad=15)
    ax.grid(True, linestyle='--', alpha=0.3)
    
    #plt.tight_layout() # Optimise l'espace autour du graphique
    plt.tight_layout()
    plt.show()
    
    return fig
# ════════════════════════════════════════════════════════════════════════════
# 11. CARACTÉRISTIQUES MÉCANIQUES — BRUTE / HOMOGÉNÉISÉE / FISSURÉE
# ════════════════════════════════════════════════════════════════════════════  

def _elem_brute(elems):
    """
    Calcule (A, y_CDG, I_propre) pour chaque élément de la section brute.
 
    Chaque élément est soit un rectangle soit un trapèze.
    elems : liste de (type, b_haut, b_bas_ou_None, hauteur, y_sommet)
    """
    A_list, y_list, I_list = [], [], []
 
    for typ, b1, b2, hi_e, y_top in elems:
        if typ == "rect":
            A = b1 * hi_e
            y = y_top + hi_e / 2.0
            I = b1 * hi_e ** 3 / 12.0
        else:  # trapèze
            A = (b1 + b2) / 2.0 * hi_e
            y_loc = hi_e * (2.0 * b2 + b1) / (3.0 * (b1 + b2))
            y = y_top + y_loc
            I = (hi_e ** 3 / 36.0) * (b1**2 + 4.0*b1*b2 + b2**2) / (b1 + b2)
 
        A_list.append(A);  y_list.append(y);  I_list.append(I)
 
    return A_list, y_list, I_list
 
 
def _inertie_beton_fissure(zones, hc):
    """
    Intègre l'aire et le moment statique de la partie comprimée (y ≤ hc).
 
    zones : liste de (y_bas, y_haut, b_bas, b_haut) de bas en haut
    hc    : hauteur de la zone comprimée (depuis la fibre inférieure)
    """
    S_fiss  = 0.0
    MS_fiss = 0.0
 
    for y1, y2, b_bot, b_top in zones:
        y_end = min(y2, hc)
        if y_end <= y1:
            continue                    # zone entièrement tendue
 
        h_z   = y2 - y1
        h_eff = y_end - y1             # hauteur comprimée dans ce segment
 
        # Largeurs aux bornes de la partie comprimée
        b_bas  = b_bot + (b_top - b_bot) * (y1   - y1) / h_z if h_z > 0 else b_bot
        b_haut = b_bot + (b_top - b_bot) * (h_eff      ) / h_z if h_z > 0 else b_top
 
        area_z = h_eff * (b_bas + b_haut) / 2.0
        if (b_bas + b_haut) > 0:
            v_loc = (h_eff / 3.0) * (b_bas + 2.0 * b_haut) / (b_bas + b_haut)
        else:
            v_loc = h_eff / 2.0
 
        S_fiss  += area_z
        MS_fiss += area_z * (y1 + v_loc)
 
    return S_fiss, MS_fiss
 
 
@func
def inerties_section_I(b, h, bs, hs, gs, bi, hi, gi,
                        asup, ainf, esup, einf,
                        n, M, N, beta, hc):
    """
    Calcule les caractéristiques mécaniques de la section en I.
 
    Retourne un dictionnaire avec :
      Section brute        : S_B, v_B, I_B
      Section homogénéisée : S_H, v_H, I_H
      Section fissurée     : S_F, v_F, I_F
    """
    h_web = h - hs - gs - gi - hi   # hauteur de l'âme
 
    # ── Description des éléments ──────────────────────────────────────────
    #    (type, b_haut, b_bas, hauteur, y_sommet)
    elems = [
        ("rect", bs,  None, hs,    0            ),  # table supérieure
        ("trap", bs,  b,    gs,    hs           ),  # congé supérieur
        ("rect", b,   None, h_web, hs + gs      ),  # âme
        ("trap", b,   bi,   gi,    hs+gs+h_web  ),  # congé inférieur
        ("rect", bi,  None, hi,    h - hi       ),  # table inférieure
    ]
 
    # ── Section brute ─────────────────────────────────────────────────────
    A_list, y_list, I_list = _elem_brute(elems)
 
    S_B = sum(A_list)
    v_B = sum(A * y for A, y in zip(A_list, y_list)) / S_B
    I_B = sum(I + A * (y - v_B)**2
              for I, A, y in zip(I_list, A_list, y_list))
 
    # ── Section homogénéisée ──────────────────────────────────────────────
    d_sup   = esup              # ordonnée acier sup / fibre inf
    d_inf   = h - einf          # ordonnée acier inf / fibre inf
    A_s_sup = n * asup  / 1e4
    A_s_inf = n * ainf  / 1e4
 
    S_H = S_B + A_s_sup + A_s_inf
    v_H = (
        sum(A * y for A, y in zip(A_list, y_list)) +
        A_s_sup * d_sup + A_s_inf * d_inf
    ) / S_H
 
    I_H = (
        sum(I + A * (y - v_H)**2 for I, A, y in zip(I_list, A_list, y_list)) +
        A_s_sup * (d_sup - v_H)**2 +
        A_s_inf * (d_inf - v_H)**2
    )
 
    # ── Section fissurée ──────────────────────────────────────────────────
    #    zones : (y_bas, y_haut, b_à_y_bas, b_à_y_haut)  de bas en haut
    zones = [
        (0,               hs,              bs, bs),
        (hs,              hs + gs,         bs, b ),
        (hs + gs,         hs+gs+h_web,     b,  b ),
        (hs+gs+h_web,     h - hi,          b,  bi),
        (h - hi,          h,               bi, bi),
    ]
 
    S_beton, MS_beton = _inertie_beton_fissure(zones, hc)
 
    S_F = S_beton + A_s_sup + A_s_inf
    v_F = (MS_beton + A_s_sup * d_sup + A_s_inf * d_inf) / S_F if S_F > 0 else 0.0
 
    # Inertie fissurée déduite de la courbure
    I_F = abs((M + N * (v_F - v_B)) / (_ES / n) / (beta / 1000.0))
 
    # Si la section est entièrement comprimée → on utilise I homogénéisée
    if hc >= h:
        S_F, v_F, I_F = S_H, v_H, I_H
 
    return {
        "S_B": S_B, "v_B": v_B, "I_B": I_B,
        "S_H": S_H, "v_H": v_H, "I_H": I_H,
        "S_F": S_F, "v_F": v_F, "I_F": I_F,
    }
 
 
@func
def e_inerties_section_I(b, h, bs, hs, gs, bi, hi, gi,
                          asup, ainf, esup, einf,
                          n, M, N, beta, hc, resultats):
    """
    Retourne une liste de caractéristiques sélectionnées (pour Excel).
    Si M < 0, la section est retournée (symétrie).
    """
    if M > 0:
        tout = inerties_section_I(b, h, bs, hs, gs, bi, hi, gi,
                                   asup, ainf, esup, einf,
                                   n, M, N, beta, hc)
    else:
        # Section retournée : on permute tables sup/inf et aciers
        tout = inerties_section_I(b, h, bi, hi, gi, bs, hs, gs,
                                   ainf, asup, einf, esup,
                                   n, -M, N, beta, hc)
 
    return [tout[r] for r in resultats.split(',')]




