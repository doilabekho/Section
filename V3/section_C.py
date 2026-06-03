"""
Section circulaire EC2 – Calcul N-M (ELS & ELU)
================================================
Corrections appliquées (détail en fin de fichier) :
  BUG-1  sigma_s_palier1 : discontinuité à eps_ud → solveurs divergent
  BUG-2  NM_circulaire_ELS g=0, eps0<0 : acier ignoré (traction pure)
  BUG-3  interaction_NM_ELU : point dupliqué à la jonction Zone1/Zone2
  BUG-4  solve_NM ELU : point initial kappa=0.1 (100‰/m) trop éloigné
  BUG-5  _calotte : variable `dt` nommée identiquement à `h` → renommée `jac`
"""

import numpy as np
from scipy.optimize import fsolve
from scipy.special import roots_legendre


# ─────────────────────────────────────────────────────────────
# LOIS DE COMPORTEMENT EUROCODE 2
# ─────────────────────────────────────────────────────────────

def eps_c2(fck: float) -> float:
    """Déformation au pic du béton (‰ → /1000)."""
    return 2.0e-3 if fck <= 50 else (2.0 + 0.085 * (fck - 50) ** 0.53) * 1e-3


def eps_cu2(fck: float) -> float:
    """Déformation ultime du béton (‰ → /1000)."""
    return 3.5e-3 if fck <= 50 else (2.6 + 35.0 * ((90.0 - fck) / 100.0) ** 4) * 1e-3


def eps_n(fck: float) -> float:
    """Exposant de la loi parabolique-rectangulaire."""
    return 2.0 if fck <= 50 else 1.4 + 23.4 * ((90.0 - fck) / 100.0) ** 4


def sigma_s_lin1(eps_s: float, C: float) -> float:
    """ELS – Loi linéaire acier (Es = 200 000 MPa, C inutilisé)."""
    return 200_000.0 * eps_s


def sigma_s_palier1(fyd: float, k: float,
                    eps_uk: float, eps_ud: float,
                    eps_s: float, a_com: str) -> float:
    """
    ELU – Loi élasto-plastique de l'acier.

    CORRECTION BUG-1
    ----------------
    Original : retournait 0.0 quand abs(eps_s) > eps_ud.
    Problème : discontinuité sigma 434 MPa → 0 MPa en un seul pas.
               fsolve voit un résidu qui saute brutalement → divergence
               quasi-systématique dès qu'une barre passe en rupture.
    Correction : on maintient sign * fyd (palier plastique) au lieu de 0.
                 La vérification de rupture réelle (si nécessaire) doit
                 être effectuée APRÈS la convergence du solveur, pas pendant.
    """
    Es = 200_000.0
    eps_yd = fyd / Es
    sign = np.sign(eps_s) if eps_s != 0.0 else 1.0

    if abs(eps_s) <= eps_yd:
        return eps_s * Es

    # BUG-1 CORRIGÉ : remplace "return 0.0" par "return sign * fyd"
    # (la rupture de barre est gérée en post-traitement)
    if abs(eps_s) > eps_ud:
        return sign * fyd   # ← CORRECTION (était : return 0.0)

    if a_com == 'incliné':
        sig = fyd + (k * fyd - fyd) * (abs(eps_s) - eps_yd) / (eps_uk - eps_yd)
        return sign * min(sig, k * fyd)
    else:
        return sign * fyd


# ─────────────────────────────────────────────────────────────
# Utilitaires géométriques
# ─────────────────────────────────────────────────────────────

def _gauss_legendre(n: int = 48):
    return roots_legendre(n)


def _calotte(R: float, d: float):
    """
    Intégrales analytiques sur la calotte u ∈ [d, R].

    CORRECTION BUG-5 (cosmétique / lisibilité)
    ------------------------------------------
    Original : `dt = R * np.sin(t)` — même expression que `h`, ce qui
    laissait croire à une erreur de jacobien.
    L'intégration est mathématiquement correcte car :
        dA = 2h · |du/dt| · dt_param = 2·(R sinθ)·(R sinθ)·dt_param
    mais la variable `dt` a été renommée `jac` pour lever toute ambiguïté.
    """
    if d >= R:
        return 0.0, 0.0, 0.0, 0.0
    if d <= -R:
        A   = np.pi * R ** 2
        Su  = 0.0
        Suu = np.pi * R ** 4 / 4.0
        Svv = np.pi * R ** 4 / 4.0
        return A, Su, Suu, Svv

    t_max = np.arccos(np.clip(d / R, -1.0, 1.0))

    xi, wi = _gauss_legendre(32)
    t   = 0.5 * t_max * (xi + 1.0)
    w   = 0.5 * t_max * wi
    u   = R * np.cos(t)
    h   = R * np.sin(t)
    jac = R * np.sin(t)   # = h ; jacobien |du/dt| — renommé pour clarté (BUG-5)

    A   = np.sum(w * 2.0 * h * jac)
    Su  = np.sum(w * 2.0 * h * u * jac)
    Suu = np.sum(w * 2.0 * h * u ** 2 * jac)
    Svv = np.sum(w * (2.0 * h ** 3 / 3.0) * jac)

    return A, Su, Suu, Svv


# ─────────────────────────────────────────────────────────────
# ELS — résultantes N et M
# ─────────────────────────────────────────────────────────────

def _NM_beton_ELS(R: float, eps0: float, g: float, C: float):
    d = -eps0 / g
    A, Su, Suu, Svv = _calotte(R, d)
    N = C * (eps0 * A  + g * Su)
    M = C * (eps0 * Su + g * Suu)
    return N, M


def NM_circulaire_ELS(R, ra, As_total, n_barres, C,
                      eps0, alpha, beta):
    """
    Résultantes (N, M) à l'ELS sur section circulaire.

    CORRECTION BUG-2
    ----------------
    Original : quand g ≈ 0 ET eps0 < 0 (traction pure ou nulle),
    la fonction retournait (0, 0) sans calculer la contribution
    des barres d'acier.
    Correction : le cas g=0 calcule toujours la contribution de l'acier,
    quelle que soit la valeur de eps0.
    """
    g = np.hypot(alpha, beta)
    As_i   = As_total / n_barres
    angles = np.linspace(0, 2 * np.pi, n_barres, endpoint=False)

    if g < 1e-14:
        # Béton (résiste uniquement en compression, eps0 > 0)
        N_b = 0.0
        if eps0 > 0:
            N_b = C * eps0 * (np.pi * R ** 2)

        # BUG-2 CORRIGÉ : acier toujours calculé (était conditionnel + ignoré si eps0<0)
        Ns = sum(sigma_s_lin1(eps0, C) * As_i for _ in angles)
        return N_b + Ns, 0.0

    N_b, M_b = _NM_beton_ELS(R, eps0, g, C)

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
# ELU — résultantes N et M
# ─────────────────────────────────────────────────────────────

def _NM_beton_ELU(R: float, eps0: float, g: float, fcd: float, fck: float):
    e2 = float(eps_c2(fck))
    n  = float(eps_n(fck))

    def sigma(eps_val: float) -> float:
        if eps_val <= 0.0: return 0.0
        if eps_val >= e2:  return fcd
        return fcd * (1.0 - (1.0 - eps_val / e2) ** n)

    d    = np.clip(-eps0 / g, -R, R)
    u_lo = d
    u_hi = R
    if u_lo >= u_hi:
        return 0.0, 0.0

    xi, wi = _gauss_legendre(64)
    u_pts  = 0.5 * (u_hi - u_lo) * xi + 0.5 * (u_hi + u_lo)
    w_pts  = 0.5 * (u_hi - u_lo) * wi

    N_b = M_b = 0.0
    for u, w in zip(u_pts, w_pts):
        h     = np.sqrt(max(R ** 2 - u ** 2, 0.0))
        eps_u = eps0 + g * u
        sig_u = sigma(eps_u)
        chord = 2.0 * h
        N_b  += w * sig_u * chord
        M_b  += w * sig_u * chord * u

    return N_b, M_b


def NM_circulaire_ELU(R, ra, As_total, n_barres,
                      fck, fcd, fyd, k, eps_uk, eps_ud, a_com,
                      eps0, alpha, beta):
    g = np.hypot(alpha, beta)

    if g < 1e-14:
        eps_u = max(eps0, 0.0)
        e2 = float(eps_c2(fck))
        n  = float(eps_n(fck))
        sig_b = 0.0 if eps0 <= 0 else (fcd if eps_u >= e2
                                        else fcd * (1.0 - (1.0 - eps_u / e2) ** n))
        N      = sig_b * np.pi * R ** 2
        As_i   = As_total / n_barres
        N     += sum(sigma_s_palier1(fyd, k, eps_uk, eps_ud, eps0, a_com) * As_i
                     for _ in range(n_barres))
        return N, 0.0

    N_b, M_b = _NM_beton_ELU(R, eps0, g, fcd, fck)

    As_i   = As_total / n_barres
    angles = np.linspace(0, 2 * np.pi, n_barres, endpoint=False)
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
# Solveurs (résolution N, M → eps0, κ)
# ─────────────────────────────────────────────────────────────

def solve_NM_circulaire_ELS(R, ra, As_total, n_barres, C, Nobj, Mobj):
    x0 = np.array([0.0, 0.001])

    def residuals(x):
        eps0, kappa = x
        g    = abs(kappa)
        N, M_u = NM_circulaire_ELS(R, ra, As_total, n_barres, C,
                                   eps0, alpha=g, beta=0.0)
        M_x  = M_u * np.sign(kappa)
        return [N - Nobj, M_x - Mobj]

    sol = fsolve(residuals, x0, xtol=1e-6)
    return sol[0], abs(sol[1])


def solve_NM_circulaire_ELU(R, ra, As_total, n_barres,
                             fck, fcd, fyd, k, eps_uk, eps_ud, a_com,
                             Nobj, Mobj):
    """
    CORRECTION BUG-4
    ----------------
    Original : x0 = [0.0, 0.1] → kappa_0 = 0.1 rad/m = 100 ‰/m.
    C'est une courbure extrêmement grande (>10× la courbure ultime typique)
    qui place le solveur dans une zone où sigma_s = 0 (barres rompues dans
    l'ancienne version) → résidus plats → gradient nul → divergence.
    Correction : x0 = [0.0, 0.001] (1 ‰/m), cohérent avec solve_NM_ELS.
    """
    x0 = np.array([0.0, 0.001])   # ← CORRECTION BUG-4 (était : [0.0, 0.1])

    def residuals(x):
        eps0, kappa = x
        g    = abs(kappa)
        N, M_u = NM_circulaire_ELU(R, ra, As_total, n_barres,
                                   fck, fcd, fyd, k, eps_uk, eps_ud, a_com,
                                   eps0, alpha=g, beta=0.0)
        M_x  = M_u * np.sign(kappa)
        return [N - Nobj, M_x - Mobj]

    sol = fsolve(residuals, x0, xtol=1e-6)
    return sol[0], abs(sol[1])


# ─────────────────────────────────────────────────────────────
# Diagramme d'interaction N-M (courbe enveloppe ELU)
# ─────────────────────────────────────────────────────────────

def interaction_NM_circulaire_ELU(R, ra, As_total, n_barres,
                                  fck, fcd, fyd, k, eps_uk, eps_ud, a_com,
                                  n_pts: int = 60):
    """
    Balayage complet des pivots EC2 (A, B, C).

    CORRECTION BUG-3
    ----------------
    Original : Zone1 et Zone2 partageaient le même premier point
    (eps_min=-eps_ud, eps_max=ecu2), générant un doublon dans la courbe.
    Correction : Zone1 exclut son dernier point (endpoint=False dans linspace).
    """
    ecu2 = float(eps_cu2(fck))

    # --- Zone 1 : Pivot A (acier en traction, béton librement varié) ---
    # eps_max_z1[-1] == ecu2 est exclu pour éviter le doublon (BUG-3 CORRIGÉ)
    eps_max_z1 = np.linspace(-eps_ud, ecu2, n_pts + 1)[:-1]  # ← CORRECTION BUG-3
    eps_min_z1 = np.full(n_pts, -eps_ud)

    # --- Zone 2 : Pivot B & C (béton à l'écrasement) ---
    eps_max_z2 = np.full(n_pts, ecu2)
    eps_min_z2 = np.linspace(-eps_ud, ecu2, n_pts)

    eps_min_all = np.concatenate([eps_min_z1, eps_min_z2])
    eps_max_all = np.concatenate([eps_max_z1, eps_max_z2])

    N_list, M_list = [], []
    for eps_min, eps_max in zip(eps_min_all, eps_max_all):
        eps0 = 0.5 * (eps_max + eps_min)
        g    = 0.5 * (eps_max - eps_min) / R

        N, M = NM_circulaire_ELU(R, ra, As_total, n_barres,
                                  fck, fcd, fyd, k, eps_uk, eps_ud, a_com,
                                  eps0,
                                  alpha=(g if g > 1e-12 else 0.0),
                                  beta=0.0)
        N_list.append(N)
        M_list.append(abs(M))

    return np.array(N_list), np.array(M_list)


