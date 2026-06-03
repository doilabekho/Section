import numpy as np
from scipy.optimize import fsolve
from scipy.special import roots_legendre

# ─────────────────────────────────────────────────────────────
# STUBS / LOIS DE COMPORTEMENT EUROCODE 2 (Ajoutés pour conformité)
# ─────────────────────────────────────────────────────────────

from beton import *   # ← fonctionne en local et sur GitHub
from acier import *   # ← fonctionne en local et sur GitHub	

# ─────────────────────────────────────────────────────────────
# Utilitaires géométriques
# ─────────────────────────────────────────────────────────────

def _gauss_legendre(n=48):
    xi, wi = roots_legendre(n)
    return xi, wi

def _calotte(R, d):
    """
    Intégrales analytiques sur la calotte u in [d, R], v in [-h, h]
    """
    if d >= R:
        return 0.0, 0.0, 0.0, 0.0
    if d <= -R:
        A   = np.pi * R**2
        Su  = 0.0
        Suu = np.pi * R**4 / 4.0
        Svv = np.pi * R**4 / 4.0
        return A, Su, Suu, Svv

    t_max = np.arccos(np.clip(d / R, -1.0, 1.0))

    xi, wi = _gauss_legendre(32)
    t  = 0.5 * t_max * (xi + 1.0)
    w  = 0.5 * t_max * wi
    u  = R * np.cos(t)
    h  = R * np.sin(t)          
    dt = R * np.sin(t)          

    A   = np.sum(w * 2.0 * h * dt)
    Su  = np.sum(w * 2.0 * h * u * dt)
    Suu = np.sum(w * 2.0 * h * u**2 * dt)
    Svv = np.sum(w * (2.0 * h**3 / 3.0) * dt)

    return A, Su, Suu, Svv


# ─────────────────────────────────────────────────────────────
# ELS — analytique exact
# ─────────────────────────────────────────────────────────────

def _NM_beton_ELS(R, eps0, g, C):
    d = -eps0 / g
    A, Su, Suu, Svv = _calotte(R, d)

    N = C * (eps0 * A   + g * Su)
    M = C * (eps0 * Su  + g * Suu)
    return N, M


def NM_circulaire_ELS(R, ra, As_total, n_barres, C, eps0, alpha, beta):
    g = np.hypot(alpha, beta)

    if g < 1e-14:
        if eps0 < 0:
            return 0.0, 0.0
        A = np.pi * R**2
        N = C * eps0 * A
        M = 0.0
        
        # CORRECTION 1 : Utilisation de la vraie loi acier à l'ELS au lieu de 'C' béton
        As_i = As_total / n_barres
        sig_s = sigma_s_lin1(eps0, C)
        N += sig_s * As_total
        return N, M

    N_b, M_b = _NM_beton_ELS(R, eps0, g, C)

    As_i   = As_total / n_barres
    angles = np.linspace(0, 2*np.pi, n_barres, endpoint=False)
    Ns = Ms = 0.0
    for a in angles:
        xs    = ra * np.cos(a)
        ys    = ra * np.sin(a)
        u_s   = (alpha * xs + beta * ys) / g     
        eps_s = eps0 + g * u_s
        sig_s = sigma_s_lin1(eps_s, C)
        F     = sig_s * As_i
        Ns   += F
        Ms   += F * u_s                           

    return N_b + Ns, M_b + Ms


# ─────────────────────────────────────────────────────────────
# ELU — Gauss sur u, analytique sur v
# ─────────────────────────────────────────────────────────────

def _NM_beton_ELU(R, eps0, g, fcd, fck):
    e2   = float(eps_c2(fck))
    n    = float(eps_n(fck))

    def sigma(eps_val):
        if eps_val <= 0.0:  return 0.0
        if eps_val >= e2:   return fcd  # Nettoyé et prolongé à l'infini (évite la rupture numérique)
        return fcd * (1.0 - (1.0 - eps_val / e2)**n)

    d    = -eps0 / g
    d    = np.clip(d, -R, R)
    u_lo = d
    u_hi = R

    if u_lo >= u_hi:
        return 0.0, 0.0

    xi, wi = _gauss_legendre(64)
    u_pts  = 0.5 * (u_hi - u_lo) * xi + 0.5 * (u_hi + u_lo)
    w_pts  = 0.5 * (u_hi - u_lo) * wi

    N_b = M_b = 0.0
    for u, w in zip(u_pts, w_pts):
        h     = np.sqrt(max(R**2 - u**2, 0.0))
        eps_u = eps0 + g * u
        sig_u = sigma(eps_u)
        chord = 2.0 * h                 

        N_b += w * sig_u * chord
        M_b += w * sig_u * chord * u    

    return N_b, M_b


def NM_circulaire_ELU(R, ra, As_total, n_barres,
                      fck, fcd, fyd, k, eps_uk, eps_ud, a_com,
                      eps0, alpha, beta):
    g = np.hypot(alpha, beta)

    if g < 1e-14:
        eps_u = max(eps0, 0.0)
        e2    = float(eps_c2(fck))
        n     = float(eps_n(fck))
        if eps_u >= e2:
            sig = fcd
        else:
            sig = fcd * (1.0 - (1.0 - eps_u / float(e2))**n)
        if eps0 <= 0: sig = 0.0
        
        N = sig * np.pi * R**2
        As_i   = As_total / n_barres
        for i in range(n_barres):
            sig_s = sigma_s_palier1(fyd, k, eps_uk, eps_ud, eps0, a_com)
            N    += sig_s * As_i
        return N, 0.0

    N_b, M_b = _NM_beton_ELU(R, eps0, g, fcd, fck)

    As_i   = As_total / n_barres
    angles = np.linspace(0, 2*np.pi, n_barres, endpoint=False)
    Ns = Ms = 0.0
    for a in angles:
        xs    = ra * np.cos(a)
        ys    = ra * np.sin(a)
        u_s   = (alpha * xs + beta * ys) / g
        eps_s = eps0 + g * u_s
        sig_s = sigma_s_palier1(fyd, k, eps_uk, eps_ud, eps_s, a_com)
        F     = sig_s * As_i
        Ns   += F
        Ms   += F * u_s

    return N_b + Ns, M_b + Ms


# ─────────────────────────────────────────────────────────────
# Solveurs Corrigés (Continuité C0/C1 assurée)
# ─────────────────────────────────────────────────────────────

def solve_NM_circulaire_ELS(R, ra, As_total, n_barres, C, Nobj, Mobj):
    x0 = np.array([0.0, 0.001])

    def residuals(x):
        eps0, kappa = x  # CORRECTION 2 : Changement de variable (courbure signée)
        g = abs(kappa)
        N, M_u = NM_circulaire_ELS(R, ra, As_total, n_barres, C, eps0, alpha=g, beta=0.0)
        M_x = M_u * np.sign(kappa) # Rend la fonction parfaitement lisse autour de g=0
        return [N - Nobj, M_x - Mobj]

    sol = fsolve(residuals, x0, xtol=1e-6)
    return sol[0], abs(sol[1])


def solve_NM_circulaire_ELU(R, ra, As_total, n_barres,
                             fck, fcd, fyd, k, eps_uk, eps_ud, a_com,
                             Nobj, Mobj):
    x0 = np.array([0.0, 0.1])

    def residuals(x):
        eps0, kappa = x  # CORRECTION 2 : Idem pour l'ELU
        g = abs(kappa)
        N, M_u = NM_circulaire_ELU(R, ra, As_total, n_barres,
                                  fck, fcd, fyd, k, eps_uk, eps_ud, a_com,
                                  eps0, alpha=g, beta=0.0)
        M_x = M_u * np.sign(kappa)
        return [N - Nobj, M_x - Mobj]

    sol = fsolve(residuals, x0, xtol=1e-6)
    return sol[0], abs(sol[1])


# ─────────────────────────────────────────────────────────────
# Diagramme d'interaction N-M (Courbe enveloppe complète)
# ─────────────────────────────────────────────────────────────

def interaction_NM_circulaire_ELU(R, ra, As_total, n_barres,
                                  fck, fcd, fyd, k, eps_uk, eps_ud, a_com,
                                  n_pts=60):
    """
    CORRECTION 3 : Balayage complet de tous les pivots EC2 (A, B et C)
    """
    ecu2 = float(eps_cu2(fck))

    N_list = []
    M_list = []

    # --- ZONE 1 : Pivot A (Rupture Acier en traction) ---
    # La fibre d'acier la plus tendue (à u = -ra) est fixée à -eps_ud.
    # On fait varier la déformation de la fibre supérieure de -eps_ud à ecu2.
    eps_min_z1 = np.full(n_pts, -eps_ud)
    eps_max_z1 = np.linspace(-eps_ud, ecu2, n_pts)

    # --- ZONE 2 : Pivot B & C (Écrasement du Béton) ---
    # La fibre supérieure du béton (à u = R) est fixée à ecu2.
    # On fait varier la déformation inférieure de -eps_ud à ecu2 (pure compression).
    eps_max_z2 = np.full(n_pts, ecu2)
    eps_min_z2 = np.linspace(-eps_ud, ecu2, n_pts)

    # Fusion des deux zones pour une courbe complète continue
    eps_min_all = np.concatenate([eps_min_z1, eps_min_z2])
    eps_max_all = np.concatenate([eps_max_z1, eps_max_z2])

    for eps_min, eps_max in zip(eps_min_all, eps_max_all):
        eps0 = 0.5 * (eps_max + eps_min)
        g    = 0.5 * (eps_max - eps_min) / R  
        
        if g < 1e-12:
            N, M = NM_circulaire_ELU(R, ra, As_total, n_barres, fck, fcd, fyd, k, eps_uk, eps_ud, a_com, eps0, 0.0, 0.0)
        else:
            N, M = NM_circulaire_ELU(R, ra, As_total, n_barres, fck, fcd, fyd, k, eps_uk, eps_ud, a_com, eps0, alpha=g, beta=0.0)
            
        N_list.append(N)
        M_list.append(abs(M))

    return np.array(N_list), np.array(M_list)
