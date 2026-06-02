
# ════════════════════════════════════════════════════════════════════════════
# ÉTAPE 1 — ÉTAT PERMANENT QP
# ════════════════════════════════════════════════════════════════════════════
#
# Équilibre sous G+P (loi béton différée via n = n_diff).
# Les câbles participent avec leur contrainte initiale sig_p (fournie).
# Pas de loi de câble modifiée ici : σ_câble = sig_p (imposé).
#
# Formule CDS (§IX.1.2) :
#   G = ∫ f_b(ε) dAb  +  Σ f_s(ε_si) Asi  +  Σ sig_pk Apk
#
# Le solveur trouve (ε_i, α_i, β_i) tel que N=NQP, My=MYQP, Mz=MZQP.
@func
def calculer_P_iso(polygon1, evi, p_acier, s_acier, a_com, n,
                p_pre, s_pre, sig_p,
                fpd, Ep, kp, eps_ukp, eps_udp,
):
    # --------------------------------------------------
    # 0. Transformation géométrique vers le repère CG
    # --------------------------------------------------
    # transformation_repere_cg calcule Cx, Cy AUTOMATIQUEMENT
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)    

    # ── Câbles (contrainte initiale imposée) ───────────────────────────────
    acierp = np.array(acier_G(p_pre, s_pre), dtype=float)
    x_p = acierp[:, 0] - Cx
    y_p = acierp[:, 1] - Cy
    A_p = acierp[:, 2]

    sig_pk = np.asarray(
        [v for v in np.atleast_1d(sig_p) if v not in ("", None)],
        dtype=float
    )


    # IMPORTANT : à l'état QP, les câbles exercent leur force initiale.
    # On ne calcule PAS une déformation variable ; on utilise directement sig_p.
    F_p     = sig_pk * A_p
    Nsp     = F_p.sum()
    Msyp    = (F_p * y_p).sum()
    Mszp    = (F_p * x_p).sum()

    return (float(-Nsp ),
            float(-Msyp ),
            float(-Mszp))

@func
def calculer_QP(polygon1, evi, p_acier, s_acier, a_com, n,
                p_pre, s_pre, sig_p,
                fpd, Ep, kp, eps_ukp, eps_udp,
                eps0, alpha, beta):
    """
    Calcule (N, My, Mz) à l'état permanent QP.

    Béton    : loi linéaire fissuré (module différé via n = Es/Ec_diff)
    Acier passif : loi linéaire
    Câbles   : contrainte imposée sig_p (pas de variation)

    CORRECTION vs code original :
      - Les câbles sont traités avec leur contrainte initiale imposée sig_p
        et non via une loi incrémentale (qui n'a pas de sens à l'état QP).
      - sig_p est un vecteur de contraintes (MPa), positif = traction câble.
    """
    # --------------------------------------------------
    # 0. Transformation géométrique vers le repère CG
    # --------------------------------------------------
    # transformation_repere_cg calcule Cx, Cy AUTOMATIQUEMENT
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)
  

    # ── Béton ─────────────────────────────────────────────────────────────
    def f_complet(x, y):
        sig = sigma_c_n1(eps0 + alpha*x + beta*y, n)  #  fissuré
        return np.array([sig, sig*y, sig*x])

    # ---- Béton plein
    res_c = polygone_integrate(f_complet, contour_cg)

    # ---- Béton évidements
    if not evidements_cg:      # cas N = 0
        res_v = np.array([0.0, 0.0, 0.0])
    else:
        res_v = np.zeros(3)
        for trou in evidements_cg:
            if len(trou) >= 3:
                res_v += polygone_integrate(f_complet, trou)

    # Résultantes béton
    Nc, Mc_y, Mc_z = res_c - res_v

    # ── Aciers passifs ─────────────────────────────────────────────────────
    if p_acier is None or len(p_acier) == 0:
        Ns, Msy, Msz = 0.0, 0.0, 0.0
    else:
        acier_data = np.array(acier_G(p_acier, s_acier), dtype=float)

    if acier_data.size == 0:
            Ns, Msy, Msz = 0.0, 0.0, 0.0

    else:
        acier = np.array(acier_G(p_acier, s_acier), dtype=float)

        # ⚠️ translation des aciers vers le repère CG
        x_s = acier[:, 0] - Cx
        y_s = acier[:, 1] - Cy
        areas_s = acier[:, 2]

        eps_s = eps0 + alpha * x_s + beta * y_s
        sig_s = sigma_s_lin1(eps_s, a_com)

        forces_s = sig_s * areas_s
        Ns = np.sum(forces_s)
        Msy = np.sum(forces_s * y_s)
        Msz = np.sum(forces_s * x_s)



    return (float(Ns + Nc ),
            float(Msy + Mc_y ),
            float(Msz  + Mc_z  ))


@func
def solve_GG_QP(polygon1, evi, p_acier, s_acier, a_com, n,
                p_pre, s_pre, sig_p,
                fpd, Ep, kp, eps_ukp, eps_udp,
                NQP, MYQP, MZQP):
    """
    Résout l'équilibre à l'état permanent QP → (ε_i, α_i, β_i).

    CORRECTIONS vs code original :
      - x0 adapté : on initialise avec un état de compression légère
        cohérent avec un état permanent sous précontrainte.
      - Cascade fsolve → root/hybr pour robustesse.
    """
    
    targets = np.array([NQP, MYQP, MZQP], float)

    def resid(ep):
        N, My, Mz = calculer_QP(
            polygon1, evi, p_acier, s_acier, a_com, n,
            p_pre, s_pre, sig_p,
            fpd, Ep, kp, eps_ukp, eps_udp,
            ep[0], ep[1], ep[2])
        return np.array([N, My, Mz]) - targets

    # Point de départ : compression légère + gradient
    x0 = np.array([0.5, 0.0, 1e-4])

    # Essai 1 : fsolve
    try:
        x, info, ier, _ = fsolve(resid, x0, full_output=True)
        if ier == 1 and np.max(np.abs(resid(x))) < 1e-3:
            return x
    except Exception:
        pass

    # Essai 2 : root/hybr
    try:
        sol = root(resid, x0, method='hybr', tol=1e-6)
        if sol.success and np.max(np.abs(sol.fun)) < 1e-3:
            return sol.x
    except Exception:
        pass

    # Essai 3 : root/lm
    sol = root(resid, x0, method='lm')
    return sol.x



# ════════════════════════════════════════════════════════════════════════════
# ÉTAPE 3 — ÉTAT FINAL ELS  (G + P + Q)
# ════════════════════════════════════════════════════════════════════════════
#
# Équilibre sous charges totales avec loi instantanée.
# Loi câble modifiée (CDS §IX.1.2) :
#   σ_p = ρ · [f_p(ε_p^(0) + ε+) − f_p(ε_p^(0))] + f_p(ε_p^(0))
#
# où ε+ = eps0 + α·x + β·y  est la déformation totale à l'état final.
#
# CORRECTION vs code original :
#   Le code original utilise eps_inc_p = eps_final - eps_QP (incrément QP→final)
#   Le manuel dit : ε+ = eps_final TOTAL (pas un incrément depuis QP).
#   La loi câble opère sur ε_p^(0) + ε+ (total depuis état 0),
#   pas sur ε_p^(0) + Δε (incrément depuis QP).

@func
def calculer_ELS(polygon1, evi, p_acier, s_acier, a_com, n,
                 p_pre, s_pre, sig_p,
                 fpd, Ep, kp, eps_ukp, eps_udp,
                 NQP, MYQP, MZQP, roh,
                 eps0, alpha, beta):
    """
    Calcule (N, My, Mz) à l'état final ELS — méthode retour à zéro INSTANTANE.

    Béton et aciers passifs : loi instantanée sur déformation TOTALE
    (pas incrémentale depuis QP, car on part de l'état zéro).
    Câbles : loi modifiée CDS §IX.1.2.
    """
    
    # --------------------------------------------------
    # 0. Transformation géométrique vers le repère CG
    # --------------------------------------------------
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

    if contour_cg is None:
        return 0.0, 0.0, 0.0

    # ── Étape 1 : état permanent QP ───────────────────────────────────────
    res_qp   = solve_GG_QP(polygon1, evi, p_acier, s_acier, a_com, n,
                            p_pre, s_pre, sig_p, fpd, Ep, kp,
                            eps_ukp, eps_udp, NQP, MYQP, MZQP)
    eps_i, alpha_i, beta_i = float(res_qp[0]), float(res_qp[1]), float(res_qp[2])

    acierp = np.array(acier_G(p_pre, s_pre), dtype=float)
    x_p = acierp[:, 0]
    y_p = acierp[:, 1]
    A_p = acierp[:, 2]


    # ── Étape 2 : déformation câble à l'état zéro ε_p^(0) ────────────────
    eps_p0 = _eps_cable_etat_zero(eps_i, alpha_i, beta_i, x_p, y_p, sig_p, n, Ep, Cx, Cy)

    # ── Étape 3a : béton — loi instantanée sur déformation totale ─────────
    # À l'état final, le béton part de zéro (décompression faite) et
    # subit la déformation totale eps+ = eps0 + α·x + β·y.
    # ── Béton ─────────────────────────────────────────────────────────────
   
    def f_complet(x, y):
        sig = sigma_c_n1(eps0 + alpha*x + beta*y, n)  #  fissuré
        return np.array([sig, sig*y, sig*x])

    # ---- Béton plein
    res_c = polygone_integrate(f_complet, contour_cg)

    # ---- Béton évidements
    if not evidements_cg:      # cas N = 0
        res_v = np.array([0.0, 0.0, 0.0])
    else:
        res_v = np.zeros(3)
        for trou in evidements_cg:
            if len(trou) >= 3:
                res_v += polygone_integrate(f_complet, trou)

    # Résultantes béton
    Nc, Mc_y, Mc_z = res_c - res_v

    # ── Étape 3b : aciers passifs — loi instantanée sur déformation totale ─
    if p_acier is None or len(p_acier) == 0:
        Ns, Msy, Msz = 0.0, 0.0, 0.0
    else:
        acier_data = np.array(acier_G(p_acier, s_acier), dtype=float)

    if acier_data.size == 0:
            Ns, Msy, Msz = 0.0, 0.0, 0.0

    else:
        acier = np.array(acier_G(p_acier, s_acier), dtype=float)
        
        x_s = acier[:, 0] - Cx
        y_s = acier[:, 1] -  Cy
        A_s = acier[:, 2]
        eps_s_final = eps0 + alpha*x_s + beta*y_s
        sig_s       = sigma_s_lin1(eps_s_final, a_com)
        F_s  = sig_s * A_s
        Ns   = F_s.sum();  Msy = (F_s*y_s).sum();  Msz = (F_s*x_s).sum()

    # ── Étape 3c : câbles — loi modifiée CDS §IX.1.2 ──────────────────────
    roh = np.asarray(
    [v for v in np.atleast_1d(roh) if v not in ("", None)],
    dtype=float)

    x_pc = x_p - Cx
    y_pc = y_p - Cy

    eps_plus   = eps0 + alpha*x_pc + beta*y_pc   # déformation totale à l'état final

    # f_p(ε_p^(0))              : contrainte câble à l'état zéro
    sig_p0_val = sigma_s_lin1_p(Ep, eps_p0, a_com)

    # f_p(ε_p^(0) + ε+)         : contrainte câble état final
    sig_p1_val = sigma_s_lin1_p(Ep, eps_p0 + eps_plus, a_com)

    # Loi modifiée : σ_p = ρ·(f_p1 − f_p0) + f_p0
    sig_pk = roh * (sig_p1_val - sig_p0_val) + sig_p0_val

    F_p   = sig_pk * A_p
    #F_p   = roh * sig_p1_val * A_p # test
    Nsp   = F_p.sum();  Msyp = (F_p*y_pc).sum();  Mszp = (F_p*x_pc).sum()
    # Effet isostatique
    Niso, Myiso, Mziso = calculer_P_iso(polygon1, evi, p_acier, s_acier, a_com, n,
            p_pre, s_pre, sig_p,
            fpd, Ep, kp, eps_ukp, eps_udp)
    coef =1
    return (float(Ns +  Nsp + Nc + coef * Niso),
            float(Msy +  Msyp + Mc_y + coef * Myiso),
            float(Msz +   Mszp + Mc_z + coef * Mziso ))


@func
def solve_GG_ELS(polygon1, evi, p_acier, s_acier, a_com, n,
                 p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
                 NQP, MYQP, MZQP, roh,
                 Nobj, Myobj, Mzobj):
    """
    Résout (N=Nobj, My=Myobj, Mz=Mzobj) → (ε₀, α, β) — ELS.
    Cascade : fsolve → hybr → lm.
    """
    targets = np.array([Nobj, Myobj, Mzobj], float)
    x0      = np.array([0.0, 1e-6, 1e-6], float)

    def resid(ep):
        N, My, Mz = calculer_ELS(
            polygon1, evi, p_acier, s_acier, a_com, n,
            p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
            NQP, MYQP, MZQP, roh, ep[0], ep[1], ep[2])
        return np.array([N, My, Mz]) - targets

    for solver_fn in (
        lambda: fsolve(resid, x0, full_output=True)[0],
        lambda: root(resid, x0, method='hybr', tol=1e-6).x,
        lambda: root(resid, x0, method='lm').x,
    ):
        try:
            x = solver_fn()
            if np.max(np.abs(resid(x))) < 1e-3:
                return x
        except Exception:
            pass

    return root(resid, x0, method='lm').x


# ════════════════════════════════════════════════════════════════════════════
# RÉSULTATS ELS
# ════════════════════════════════════════════════════════════════════════════

@func
def resultats_GG_ELS(polygon1, evi, p_acier, s_acier, a_com, n,
                     p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
                     NQP, MYQP, MZQP, roh,
                     eps0, alpha, beta):
    """Retourne un dict de résultats ELS pour le plan (ε₀, α, β)."""
    # --------------------------------------------------
    # 0. Transformation géométrique vers le repère CG
    # --------------------------------------------------
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

    if contour_cg is None:
        return {}

    # --------------------------------------------------
    # 1. Vérification contraintes béton (repère CG)
    # --------------------------------------------------
    if len(contour_cg) > 0:
        poly_array = np.asarray(contour_cg, dtype=float)

        eps_c = eps0 + alpha * poly_array[:, 0] + beta * poly_array[:, 1]
        espcmax = float(np.max(eps_c))
        espcmin = float(np.min(eps_c))

        sig_c_max = float(sigma_c_n1(espcmax, n))
        sig_c_min = float(sigma_c_n1(espcmin, n))
    else:
        espcmax = espcmin = 0.0
        sig_c_max = sig_c_min = 0.0

    # --------------------------------------------------
    # 2. Vérification contraintes acier (repère CG)
    # --------------------------------------------------
    if p_acier is None or len(p_acier) == 0:
        espsmax = espsmin = ""
        sig_s_max = sig_s_min = ""
    else:
        acier = acier_G(p_acier, s_acier)

        if len(acier) == 0:
            espsmax = espsmin = ""
            sig_s_max = sig_s_min = ""
        else:
            acier_array = np.asarray(acier, dtype=float)

            # coordonnées déjà TRANSLATÉES vers CG
            x_s = acier_array[:, 0] - Cx
            y_s = acier_array[:, 1] - Cy

            eps_s = eps0 + alpha * x_s + beta * y_s
            espsmax = float(np.max(eps_s))
            espsmin = float(np.min(eps_s))

            sig_s_max = float(sigma_s_lin1(espsmax, a_com))
            sig_s_min = float(sigma_s_lin1(espsmin, a_com))

    # Câbles
    res_qp = solve_GG_QP(polygon1, evi, p_acier, s_acier, a_com, n,
                          p_pre, s_pre, sig_p, fpd, Ep, kp,
                          eps_ukp, eps_udp, NQP, MYQP, MZQP)
    eps_i, alpha_i, beta_i = float(res_qp[0]), float(res_qp[1]), float(res_qp[2])

    acierp = np.array(acier_G(p_pre, s_pre), dtype=float)
    x_p = acierp[:, 0]
    y_p = acierp[:, 1]
    A_p = acierp[:, 2]

    eps_p0    = _eps_cable_etat_zero(eps_i, alpha_i, beta_i, x_p, y_p, sig_p, n, Ep, Cx, Cy)

    x_pc = x_p - Cx
    y_pc = y_p - Cy
    eps_plus  = eps0 + alpha*x_pc + beta*y_pc
    sig_p0_v  = sigma_s_lin1_p(Ep, eps_p0, a_com)
    sig_p1_v  = sigma_s_lin1_p(Ep, eps_p0 + eps_plus, a_com)
    contrai = sig_p1_v
    #contrai   = roh*(sig_p1_v - sig_p0_v) + sig_p0_v
    espp      = eps_s_palier_p(fpd, Ep, kp, eps_ukp, eps_udp, contrai)

    A_com = S_com(polygon1, evi, eps0, alpha, beta)
    N, My, Mz = calculer_ELS(
        polygon1, evi, p_acier, s_acier, a_com, n,
        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
        NQP, MYQP, MZQP, roh, eps0, alpha, beta)

    return {
        "ACOM":      A_com,
        "EPS_C_MAX": espcmax,
        "EPS_C_MIN": espcmin,
        "SIG_C_MAX": sig_c_max,
        "SIG_C_MIN": sig_c_min,
        "EPS_P_MAX": float(espp.max()),  "EPS_P_MIN": float(espp.min()),
        "EPS_S_MAX": espsmax,
        "EPS_S_MIN": espsmin,
        "SIG_S_MAX": sig_s_max,
        "SIG_S_MIN": sig_s_min,
        "SIG_P_MAX": float(contrai.max()), "SIG_P_MIN": float(contrai.min()),
        "SIG_P":     contrai.tolist(),     "EPS_P":     espp.tolist(),
        "N": float(N), "MY": float(My), "MZ": float(Mz),
    }


@func
def e_resultats_GG_ELS(polygon1, evi, p_acier, s_acier, a_com, n,
                        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
                        NQP, MYQP, MZQP, roh,
                        eps0, alpha, beta, resultats):
    tout = resultats_GG_ELS(polygon1, evi, p_acier, s_acier, a_com, n,
                            p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
                            NQP, MYQP, MZQP, roh, eps0, alpha, beta)
    return [tout[r] for r in resultats.split(',')]


# ════════════════════════════════════════════════════════════════════════════
# ÉTAPE 3 — ÉTAT FINAL ELU  (G + P + Q)
# ════════════════════════════════════════════════════════════════════════════
#
# Identique à ELS mais avec :
#   Béton    : loi parabole-rectangle (sigma_c_pararect1)
#   Acier passif : loi palier (sigma_s_palier1)
#   Câbles   : même loi modifiée CDS §IX.1.2, avec f_p = loi palier

@func
def calculer_N_My_Mz_ELU_pararect(
        polygon1, evi, p_acier, s_acier, a_com,
        fck, fcd, fyd, k, eps_uk, eps_ud,
        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
        n, NQP, MYQP, MZQP, roh,
        eps0, alpha, beta):
    """
    Calcule (N, My, Mz) ELU — méthode retour à zéro INSTANTANE.

    CORRECTIONS vs code original :
      1. ε+ est la déformation TOTALE (pas un incrément depuis QP).
      2. La loi câble utilise ε_p^(0) + ε+ (pas ε_p^(0) + Δε_QP→final).
      3. σ_b à l'état QP est corrigé : signe − dans la formule ε_p^(0).
    """
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

    if contour_cg is None:
        return 0.0, 0.0, 0.0
    # ── État QP ──────────────────────────────────────────────────────────
    res_qp   = solve_GG_QP(polygon1, evi, p_acier, s_acier, a_com, n,
                            p_pre, s_pre, sig_p, fpd, Ep, kp,
                            eps_ukp, eps_udp, NQP, MYQP, MZQP)
    eps_i, alpha_i, beta_i = float(res_qp[0]), float(res_qp[1]), float(res_qp[2])

    acierp = np.array(acier_G(p_pre, s_pre), dtype=float)
    x_p = acierp[:, 0]
    y_p = acierp[:, 1]
    A_p = acierp[:, 2]

    # ── ε_p^(0) ──────────────────────────────────────────────────────────
    eps_p0 = _eps_cable_etat_zero(eps_i, alpha_i, beta_i, x_p, y_p, sig_p, n, Ep, Cx, Cy)

    # ── Béton ELU — loi parabole-rectangle ───────────────────────────────
    def f_complet(x, y):
        sig = sigma_c_pararect1(fck, fcd, eps0 + alpha*x + beta*y)
        return np.array([sig, sig*y, sig*x])

    # ---- Béton plein
    res_c = polygone_integrate(f_complet, contour_cg)

    # ---- Béton évidements
    if not evidements_cg:      # cas N = 0
        res_v = np.array([0.0, 0.0, 0.0])
    else:
        res_v = np.zeros(3)
        for trou in evidements_cg:
            if len(trou) >= 3:
                res_v += polygone_integrate(f_complet, trou)

    # Résultantes béton
    Nc, Mc_y, Mc_z = res_c - res_v

    # ── Aciers passifs ELU — loi palier ──────────────────────────────────
    if p_acier is None or len(p_acier) == 0:
        Ns, Msy, Msz = 0.0, 0.0, 0.0
    else:
        acier_data = np.array(acier_G(p_acier, s_acier), dtype=float)

    if acier_data.size == 0:
            Ns, Msy, Msz = 0.0, 0.0, 0.0

    else:
        acier = np.array(acier_G(p_acier, s_acier), dtype=float)
        x_s = acier[:, 0] - Cx
        y_s = acier[:, 1] - Cy
        A_s = acier[:, 2]
        eps_s_final = eps0 + alpha*x_s + beta*y_s
        sig_s       = sigma_s_palier1(fyd, k, eps_uk, eps_ud, eps_s_final, a_com)
        F_s  = sig_s * A_s
        Ns   = F_s.sum();  Msy = (F_s*y_s).sum();  Msz = (F_s*x_s).sum()

    # ── Câbles ELU — loi modifiée CDS §IX.1.2 ────────────────────────────
    roh = np.asarray(
    [v for v in np.atleast_1d(roh) if v not in ("", None)],
    dtype=float)
    x_pc = x_p - Cx
    y_pc = y_p - Cy
    eps_plus   = eps0 + alpha*x_pc + beta*y_pc

    sig_p0_val = sigma_s_palier_p1(fpd, Ep, kp, eps_ukp, eps_udp, eps_p0, a_com)
    sig_p1_val = sigma_s_palier_p1(fpd, Ep, kp, eps_ukp, eps_udp, eps_p0 + eps_plus, a_com)
    sig_pk     = roh*(sig_p1_val - sig_p0_val) + sig_p0_val

    F_p   = sig_pk * A_p
    Nsp   = F_p.sum();  Msyp = (F_p*y_pc).sum();  Mszp = (F_p*x_pc).sum()

    Niso, Myiso, Mziso = calculer_P_iso(polygon1, evi, p_acier, s_acier, a_com, n,
            p_pre, s_pre, sig_p,
            fpd, Ep, kp, eps_ukp, eps_udp)
    coef = 1
    return (float(Ns +  Nsp + Nc + coef * Niso),
            float(Msy +  Msyp + Mc_y + coef * Myiso),
            float(Msz +   Mszp + Mc_z + coef * Mziso ))
@func
def solve_GG_ELU_pararect(
        polygon1, evi, p_acier, s_acier, a_com,
        fck, fcd, fyd, k, eps_uk, eps_ud,
        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
        n, NQP, MYQP, MZQP, roh,
        Nobj, Myobj, Mzobj):
    """
    Résout (N=Nobj, My=Myobj, Mz=Mzobj) → (ε₀, α, β) — ELU.
    """
    targets = np.array([Nobj, Myobj, Mzobj], float)
    x0      = np.array([0.5, 0.0, 1e-4], float)

    def resid(ep):
        N, My, Mz = calculer_N_My_Mz_ELU_pararect(
            polygon1, evi, p_acier, s_acier, a_com,
            fck, fcd, fyd, k, eps_uk, eps_ud,
            p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
            n, NQP, MYQP, MZQP, roh, ep[0], ep[1], ep[2])
        return np.array([N, My, Mz]) - targets

    for solver_fn in (
        lambda: fsolve(resid, x0, full_output=True)[0],
        lambda: root(resid, x0, method='hybr', tol=1e-6).x,
        lambda: root(resid, x0, method='lm').x,
    ):
        try:
            x = solver_fn()
            if np.max(np.abs(resid(x))) < 1e-3:
                return x
        except Exception:
            pass

    return root(resid, x0, method='lm').x


# ════════════════════════════════════════════════════════════════════════════
# RÉSULTATS ELU
# ════════════════════════════════════════════════════════════════════════════

@func
def resultats_GG_ELU_pararect(
        polygon1, evi, p_acier, s_acier, a_com,
        fck, fcd, fyd, k, eps_uk, eps_ud,
        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
        n, NQP, MYQP, MZQP, roh,
        eps0, alpha, beta):
    """Retourne un dict de résultats ELU pour le plan (ε₀, α, β)."""
    contour_cg, evidements_cg, Cx, Cy = transformation_repere_cg(polygon1, evi)

     # --------------------------------------------------
    # 1. Vérification contraintes béton (repère CG)
    # --------------------------------------------------
    if len(contour_cg) > 0:
        poly_array = np.asarray(contour_cg, dtype=float)

        eps_c = eps0 + alpha * poly_array[:, 0] + beta * poly_array[:, 1]
        espcmax = float(np.max(eps_c))
        espcmin = float(np.min(eps_c))

        sig_c_max = float(sigma_c_pararect1(fck, fcd, espcmax))
        sig_c_min = float(sigma_c_pararect1(fck, fcd, espcmin))
    else:
        espcmax = espcmin = 0.0
        sig_c_max = sig_c_min = 0.0

    # --------------------------------------------------
    # 2. Vérification contraintes acier (repère CG)
    # --------------------------------------------------
    if p_acier is None or len(p_acier) == 0:
        epssmax = epssmin = ""
        sig_s_max = sig_s_min = ""
        pour_a_pl = 0.0
    else:
        acier = acier_G(p_acier, s_acier)

        if len(acier) == 0:
            epssmax = epssmin = ""
            sig_s_max = sig_s_min = ""
            pour_a_pl = 0.0
        else:
            acier_array = np.asarray(acier, dtype=float)

            # coordonnées déjà TRANSLATÉES vers CG
            x_s = acier_array[:, 0] - Cx
            y_s = acier_array[:, 1] - Cy

            eps_s = eps0 + alpha * x_s + beta * y_s
            epssmax = float(np.max(eps_s))
            epssmin = float(np.min(eps_s))

            # Pourcentage plastifié (Es ≈ 200 000 MPa)
            plast = np.sum(np.abs(eps_s) >= (float(fyd) / 200.0))
            pour_a_pl = float(plast / len(eps_s) * 100.0)

            sig_s_max = float(sigma_s_palier1(
                fyd, k, eps_uk, eps_ud, epssmax, a_com
            ))
            sig_s_min = float(sigma_s_palier1(
                fyd, k, eps_uk, eps_ud, epssmin, a_com
            ))

    # Câbles
    res_qp = solve_GG_QP(polygon1, evi, p_acier, s_acier, a_com, n,
                          p_pre, s_pre, sig_p, fpd, Ep, kp,
                          eps_ukp, eps_udp, NQP, MYQP, MZQP)
    eps_i, alpha_i, beta_i = float(res_qp[0]), float(res_qp[1]), float(res_qp[2])

    acierp = np.array(acier_G(p_pre, s_pre), dtype=float)
    x_p = acierp[:, 0]
    y_p = acierp[:, 1]
    A_p = acierp[:, 2]

    eps_p0    = _eps_cable_etat_zero(eps_i, alpha_i, beta_i, x_p, y_p, sig_p, n, Ep, Cx, Cy)

    x_pc = x_p - Cx
    y_pc = y_p - Cy
    eps_plus = eps0 + alpha*x_pc + beta*y_pc
    sig_p0_v = sigma_s_palier_p1(fpd, Ep, kp, eps_ukp, eps_udp, eps_p0, a_com)
    sig_p1_v = sigma_s_palier_p1(fpd, Ep, kp, eps_ukp, eps_udp, eps_p0 + eps_plus, a_com)
    contrai  = sig_p1_v
    #contrai  = roh*(sig_p1_v - sig_p0_v) + sig_p0_v
    espp     = eps_s_palier_p(fpd, Ep, kp, eps_ukp, eps_udp, contrai)

    # Plastification acier passif


    A_com = S_com(polygon1, evi, eps0, alpha, beta)
    N, My, Mz = calculer_N_My_Mz_ELU_pararect(
        polygon1, evi, p_acier, s_acier, a_com,
        fck, fcd, fyd, k, eps_uk, eps_ud,
        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
        n, NQP, MYQP, MZQP, roh, eps0, alpha, beta)

    return {
        "EPS0": eps0, "ALPHA": alpha, "BETA": beta, "ACOM": A_com,
        "EPS_C_MAX": espcmax,
        "EPS_C_MIN": espcmin,
        "SIG_C_MAX": sig_c_max,
        "SIG_C_MIN": sig_c_min,
        "EPS_P_MAX": float(espp.max()),  "EPS_P_MIN": float(espp.min()),
        "EPS_S_MAX": epssmax,
        "EPS_S_MIN": epssmin,
        "SIG_S_MAX": sig_s_max,
        "SIG_S_MIN": sig_s_min,
        "SIG_P_MAX": float(contrai.max()), "SIG_P_MIN": float(contrai.min()),
        "SIG_P":     contrai.tolist(),     "EPS_P":     espp.tolist(),
        "N": float(N), "MY": float(My), "MZ": float(Mz),
        "PA": pour_a_pl,
    }


@func
def e_resultats_GG_ELU_pararect(
        polygon1, evi, p_acier, s_acier, a_com,
        fck, fcd, fyd, k, eps_uk, eps_ud,
        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
        n, NQP, MYQP, MZQP, roh,
        eps0, alpha, beta, resultats):
    tout = resultats_GG_ELU_pararect(
        polygon1, evi, p_acier, s_acier, a_com,
        fck, fcd, fyd, k, eps_uk, eps_ud,
        p_pre, s_pre, sig_p, fpd, Ep, kp, eps_ukp, eps_udp,
        n, NQP, MYQP, MZQP, roh, eps0, alpha, beta)
    return [tout[r] for r in resultats.split(',')]
