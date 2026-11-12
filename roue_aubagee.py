#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Roue aubagée axiale simplifiée (disque + pales monobloc), maillage 100 % HEXA8.

Principe
--------
1. Gmsh maille en quadrangles l'empreinte d'UN secteur sur la jante, dans le
   plan développé (s = R_jante * theta, z). Le secteur est un parallélogramme
   incliné selon le calage de la pale ; il contient le profil cambré à bords
   arrondis (surface "pale") et le reste de la jante (surface "hors pale").
2. Cette empreinte 2D est balayée radialement :
     - vers l'intérieur -> disque (alésage / voile / jante, largeur axiale w(r))
     - vers l'extérieur -> pale (uniquement les quads du profil, avec vrillage)
   Les noeuds sont partagés à l'emplanture : maillage conforme pale/disque.
3. Le secteur est recopié N fois par rotation, les noeuds des faces de
   périodicité sont fusionnés -> roue complète 360° conforme.
4. Le maillage est injecté dans Gmsh (entités discrètes + groupes physiques)
   et écrit au format MED.

Axe de rotation : Z. Unités : mm. Amont = z négatif (bord d'attaque).

Groupes MED
-----------
  Volumes  : DISQUE, PALES
  Surfaces : ALESAGE, AMONT, AVAL, JANTE, PEAU_PALES, TETE_PALES
             (quadrangles orientés normale sortante)
  Point    : CENTRE (maille POI1 sur l'axe, en z = Z_CENTRE) : noeud maître
             à 6 ddl (élément discret) à relier à ALESAGE dans le calcul
             pour encastrer / piloter la roue.

Usage : python3 roue_aubagee.py [--gui]
"""
import sys
import numpy as np
import gmsh

# =============================================================================
# PARAMÈTRES
# =============================================================================
# --- Roue ---
N_PALES = 20            # nombre de pales (= nombre de secteurs)

# --- Disque (profil méridien : alésage / voile / jante) ---
R_ALESAGE = 40.0        # rayon intérieur
R_MOYEU = 60.0          # fin du moyeu (largeur constante W_MOYEU)
R_VOILE_INT = 75.0      # début du voile (largeur W_VOILE)
R_VOILE_EXT = 125.0     # fin du voile
R_JANTE_INT = 138.0     # début de la jante (largeur W_JANTE)
R_JANTE = 150.0         # rayon extérieur de la jante = pied de pale
W_MOYEU = 40.0          # largeurs axiales
W_VOILE = 16.0
W_JANTE = 46.0

# --- Pale ---
HAUTEUR = 60.0          # hauteur de pale (envergure)
CORDE = 40.0            # corde en pied
EPAISSEUR = 4.0         # épaisseur max du profil
R_BA = 1.2              # rayon du bord d'attaque (demi-cercle)
R_BF = 0.9              # rayon du bord de fuite (demi-cercle)
CAMBRURE = 0.08         # flèche de cambrure / corde
POS_EPAISSEUR = 0.40    # position de l'épaisseur max (fraction de corde)
CALAGE_DEG = 30.0       # calage en pied (angle corde / axe Z)
VRILLAGE_DEG = 20.0     # calage supplémentaire en tête (linéaire en envergure)
EFFILEMENT = 0.85       # corde tête / corde pied

# --- Tailles de maille (tailles finales visées) ---
LC_PALE = 2.0           # le long de la corde et au voisinage du profil
N_EP_COEUR = 4          # éléments dans l'épaisseur du coeur du profil (pair)
N_EP_COURONNE = 2       # éléments dans la couronne, de chaque côté (pair)
FRACTION_COURONNE = 0.50  # part de la demi-épaisseur occupée par la couronne
EP_TRANSITION = 2.5     # épaisseur de la couronne structurée autour du profil
N_TRANSITION = 2        # éléments dans cette couronne (pair)
LC_JANTE = 6.0          # loin de la pale, sur la jante
LC_RADIAL_DISQUE = 6.0  # taille radiale dans le disque
LC_RADIAL_PIED = 1.5    # taille radiale juste sous le pied de pale
RAISON_RADIALE = 1.3    # progression géométrique entre les deux
N_COUCHES_PALE = 24     # nombre de couches d'éléments sur l'envergure
N_LISSAGE = 50          # itérations de lissage de la zone non structurée
DETENTE_ITER = 1000      # détente de l'empreinte dans le disque (0 = désactivée)

# --- Noeud maître sur l'axe ---
Z_CENTRE = 0.0          # position axiale du noeud CENTRE (0 = milieu du disque)

FICHIER = "roue_aubagee.med"


# =============================================================================
# 1. Empreinte 2D d'un secteur (plan développé s, z) maillée en quads par Gmsh
# =============================================================================
def _vers_sz(X, Y):
    """Repère profil (X le long de la corde, Y normal, origine à mi-corde)
    -> plan développé (s, z). Corde dirigée selon (sin g, cos g), BA côté z < 0."""
    g = np.radians(CALAGE_DEG)
    s = CORDE * (X * np.sin(g) + Y * np.cos(g))
    z = CORDE * (X * np.cos(g) - Y * np.sin(g))
    return np.column_stack([np.atleast_1d(s), np.atleast_1d(z)])


def _ligne_moyenne(x):
    """Ligne moyenne parabolique : point, tangente et normale unitaires."""
    x = np.asarray(x, dtype=float)
    yc = 4.0 * CAMBRURE * x * (1.0 - x)
    dyc = 4.0 * CAMBRURE * (1.0 - 2.0 * x)
    nrm = np.hypot(1.0, dyc)
    return x, yc, 1.0 / nrm, dyc / nrm          # x, y, tx, ty  (n = (-ty, tx))


# stations des centres des bords arrondis (en fraction de corde)
X_BA = R_BA / CORDE
X_BF = 1.0 - R_BF / CORDE


def demi_epaisseur(x):
    """Demi-épaisseur / corde entre les deux centres de bord : R_BA au bord
    d'attaque, EPAISSEUR/2 en POS_EPAISSEUR, R_BF au bord de fuite, avec
    pentes nulles aux extrémités (raccord tangent avec les demi-cercles)."""
    x = np.asarray(x, dtype=float)
    hm, ha, hb = 0.5 * EPAISSEUR / CORDE, R_BA / CORDE, R_BF / CORDE
    xm = POS_EPAISSEUR
    t1 = np.clip((x - X_BA) / (xm - X_BA), 0.0, 1.0)
    t2 = np.clip((x - xm) / (X_BF - xm), 0.0, 1.0)
    av = ha + (hm - ha) * 0.5 * (1.0 - np.cos(np.pi * t1))
    ar = hm + (hb - hm) * 0.5 * (1.0 - np.cos(np.pi * t2))
    return np.where(x <= xm, av, ar)


def point_profil(x, eta, dec=0.0):
    """Point du corps du profil : station x (X_BA..X_BF), eta = position dans
    l'épaisseur (-1 intrados ... +1 extrados). dec = décalage normal vers
    l'extérieur (mm), pour la couronne de transition côté jante."""
    x, yc, tx, ty = _ligne_moyenne(x)
    h = demi_epaisseur(x) + dec / CORDE
    return _vers_sz(x - eta * h * ty - 0.5, yc + eta * h * tx)


def point_bord(bord, angle, dec=0.0):
    """Point du demi-cercle de bord ('ba' ou 'bf'). angle (rad) mesuré depuis
    la pointe du bord : +pi/2 = raccord extrados, -pi/2 = raccord intrados."""
    x0, r, sens = (X_BA, R_BA, -1.0) if bord == "ba" else (X_BF, R_BF, +1.0)
    x, yc, tx, ty = _ligne_moyenne(x0)
    r = (r + dec) / CORDE
    ct, st = sens * r * np.cos(angle), r * np.sin(angle)
    return _vers_sz(x + ct * tx - st * ty - 0.5, yc + ct * ty + st * tx)


def empreinte_2d():
    """Renvoie (noeuds (n,2), quads (m,4) antihoraires, masque 'quad de pale').

    Le profil est maillé en blocs structurés (topologie en O) :
    un coeur + une couronne (extrados, intrados, calotte BA, calotte BF)."""
    pas = 2.0 * np.pi * R_JANTE / N_PALES
    tg = np.tan(np.radians(CALAGE_DEG))
    zmin, zmax = -0.5 * W_JANTE, 0.5 * W_JANTE

    # Le contour est découpé en 4 courbes ; les points de découpe sont sur les
    # demi-cercles de bord, à +/- ANG de la pointe, pour que les rayons de la
    # couronne partent en diagonale vers les coins du coeur.
    ANG = np.radians(65.0)
    eta_c = 1.0 - FRACTION_COURONNE                # frontière du coeur
    xs = np.linspace(X_BA, X_BF, 61)
    arc = np.linspace(ANG, 0.5 * np.pi, 7)
    cal = np.linspace(ANG, -ANG, 25)
    contour = {
        "ext": np.vstack([point_bord("ba", arc)[:-1], point_profil(xs, +1.0),
                          point_bord("bf", arc[::-1])[1:]]),
        "int": np.vstack([point_bord("ba", -arc)[:-1], point_profil(xs, -1.0),
                          point_bord("bf", -arc[::-1])[1:]]),
        "ba": point_bord("ba", cal),                # extrados -> intrados
        "bf": point_bord("bf", cal),
    }
    xc = xs
    d = EP_TRANSITION                              # contour décalé (transition)
    decale = {
        "ext": np.vstack([point_bord("ba", arc, d)[:-1], point_profil(xs, +1.0, d),
                          point_bord("bf", arc[::-1], d)[1:]]),
        "int": np.vstack([point_bord("ba", -arc, d)[:-1], point_profil(xs, -1.0, d),
                          point_bord("bf", -arc[::-1], d)[1:]]),
        "ba": point_bord("ba", cal, d),
        "bf": point_bord("bf", cal, d),
    }
    tout = np.vstack(list(decale.values()))
    sp, zp = tout[:, 0], tout[:, 1]
    # garde-fous géométriques
    if zp.min() < zmin + LC_PALE or zp.max() > zmax - LC_PALE:
        raise ValueError("Le profil dépasse de la jante : augmenter W_JANTE.")
    marge = 0.5 * pas - np.abs(sp - zp * tg).max()
    if marge < 1.5 * LC_PALE:
        raise ValueError("Le profil ne tient pas dans le secteur : réduire "
                         "N_PALES, la corde ou l'épaisseur.")

    gmsh.model.add("empreinte")
    occ = gmsh.model.occ
    # la subdivision finale divise la taille par 2 -> on demande le double
    f = 2.0

    # parallélogramme du secteur (côtés inclinés du calage => périodiques)
    P = [occ.addPoint(-0.5 * pas + zmin * tg, zmin, 0),
         occ.addPoint(+0.5 * pas + zmin * tg, zmin, 0),
         occ.addPoint(+0.5 * pas + zmax * tg, zmax, 0),
         occ.addPoint(-0.5 * pas + zmax * tg, zmax, 0)]
    L = [occ.addLine(P[i], P[(i + 1) % 4]) for i in range(4)]

    def spline(xy, p_deb, p_fin):
        mil = [occ.addPoint(a, b, 0) for a, b in xy[1:-1]]
        return occ.addSpline([p_deb] + mil + [p_fin])

    def pt(xy):
        return occ.addPoint(xy[0], xy[1], 0)

    # sommets : contour (majuscules) et coeur (minuscules)
    Au, Bu = pt(contour["ext"][0]), pt(contour["ext"][-1])
    Al, Bl = pt(contour["int"][0]), pt(contour["int"][-1])
    coeur_u, coeur_l = point_profil(xc, +eta_c), point_profil(xc, -eta_c)
    au, bu, al, bl = pt(coeur_u[0]), pt(coeur_u[-1]), pt(coeur_l[0]), pt(coeur_l[-1])
    # courbes du contour
    c_ext = spline(contour["ext"], Au, Bu)
    c_int = spline(contour["int"], Al, Bl)
    c_ba = spline(contour["ba"], Au, Al)
    c_bf = spline(contour["bf"], Bu, Bl)
    # courbes du coeur
    k_u, k_l = spline(coeur_u, au, bu), spline(coeur_l, al, bl)
    k_ba, k_bf = occ.addLine(au, al), occ.addLine(bu, bl)
    # rayons de la couronne
    rAu, rAl = occ.addLine(Au, au), occ.addLine(Al, al)
    rBu, rBl = occ.addLine(Bu, bu), occ.addLine(Bl, bl)

    def surf(courbes):
        return occ.addPlaneSurface([occ.addCurveLoop(courbes)])

    s_blocs = [surf([k_u, k_bf, k_l, k_ba]),       # coeur
               surf([c_ext, rBu, k_u, rAu]),       # couronne extrados
               surf([c_int, rBl, k_l, rAl]),       # couronne intrados
               surf([c_ba, rAl, k_ba, rAu]),       # calotte BA
               surf([c_bf, rBl, k_bf, rBu])]       # calotte BF
    # couronne de transition structurée autour du profil (côté jante)
    Du, Eu = pt(decale["ext"][0]), pt(decale["ext"][-1])
    Dl, El = pt(decale["int"][0]), pt(decale["int"][-1])
    d_ext, d_int = spline(decale["ext"], Du, Eu), spline(decale["int"], Dl, El)
    d_ba, d_bf = spline(decale["ba"], Du, Dl), spline(decale["bf"], Eu, El)
    tAu, tAl = occ.addLine(Au, Du), occ.addLine(Al, Dl)
    tBu, tBl = occ.addLine(Bu, Eu), occ.addLine(Bl, El)
    s_trans = [surf([c_ext, tBu, d_ext, tAu]), surf([c_int, tBl, d_int, tAl]),
               surf([c_ba, tAl, d_ba, tAu]), surf([c_bf, tBl, d_bf, tBu])]
    b_ext = occ.addCurveLoop(L)
    b_prof = occ.addCurveLoop([d_ext, d_bf, d_int, d_ba])
    s_libre = occ.addPlaneSurface([b_ext, b_prof])
    occ.synchronize()

    # discrétisation structurée du profil (nombres finaux pairs, cf. subdivision)
    long_ext = gmsh.model.occ.getMass(1, c_ext)
    n_corde = max(2, int(round(long_ext / LC_PALE / 2.0)))      # demi-nombre
    n_coeur = max(1, N_EP_COEUR // 2)
    n_cour = max(1, N_EP_COURONNE // 2)
    n_trans = max(1, N_TRANSITION // 2)
    for c in (tAu, tAl, tBu, tBl):
        gmsh.model.mesh.setTransfiniteCurve(c, n_trans + 1)
    for c in (c_ext, c_int, k_u, k_l, d_ext, d_int):
        gmsh.model.mesh.setTransfiniteCurve(c, n_corde + 1)
    for c in (c_ba, c_bf, k_ba, k_bf, d_ba, d_bf):
        gmsh.model.mesh.setTransfiniteCurve(c, n_coeur + 1)
    for c in (rAu, rAl, rBu, rBl):
        gmsh.model.mesh.setTransfiniteCurve(c, n_cour + 1)
    for s in s_blocs + s_trans:
        gmsh.model.mesh.setTransfiniteSurface(s)
    c_prof = [d_ext, d_bf, d_int, d_ba]

    # même discrétisation sur les deux côtés périodiques
    n_cote = max(3, int(round(W_JANTE / np.cos(np.radians(CALAGE_DEG))
                              / (f * LC_JANTE))) + 1)
    for c in (L[1], L[3]):
        gmsh.model.mesh.setTransfiniteCurve(c, n_cote)

    # champ de taille : fin près du profil, grossier loin
    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "CurvesList", c_prof)
    gmsh.model.mesh.field.setNumber(fd, "Sampling", 400)
    ft = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(ft, "InField", fd)
    gmsh.model.mesh.field.setNumber(ft, "SizeMin", f * LC_PALE)
    gmsh.model.mesh.field.setNumber(ft, "SizeMax", f * LC_JANTE)
    gmsh.model.mesh.field.setNumber(ft, "DistMin", 0.5 * LC_PALE)
    gmsh.model.mesh.field.setNumber(ft, "DistMax", 0.3 * pas)
    gmsh.model.mesh.field.setAsBackgroundMesh(ft)
    for opt in ("MeshSizeFromPoints", "MeshSizeFromCurvature",
                "MeshSizeExtendFromBoundary"):
        gmsh.option.setNumber("Mesh." + opt, 0)

    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)     # blossom
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)   # garantit 100 % quads
    gmsh.option.setNumber("Mesh.Smoothing", 20)

    def quads(surf, renum):
        types, _, _ = gmsh.model.mesh.getElements(2, surf)
        if list(types) != [3]:
            raise RuntimeError("Le maillage 2D n'est pas 100 % quadrangles.")
        _, nt = gmsh.model.mesh.getElementsByType(3, surf)
        return renum[nt.astype(np.int64)].reshape(-1, 4)

    def extraire():
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        xyz = xyz.reshape(-1, 3)[:, :2]
        renum = np.full(int(tags.max()) + 1, -1, dtype=np.int64)
        renum[tags.astype(np.int64)] = np.arange(len(tags))
        q_pale = np.vstack([quads(s, renum) for s in s_blocs])
        q_trans = np.vstack([quads(s, renum) for s in s_trans])
        q = np.vstack([q_pale, q_trans, quads(s_libre, renum)])
        n_struct = len(q_pale) + len(q_trans)
        est_pale = np.zeros(len(q), dtype=bool)
        est_pale[:len(q_pale)] = True

        # ne garder que les noeuds utilisés, orienter les quads en antihoraire
        util = np.unique(q)
        renum2 = np.full(len(xyz), -1, dtype=np.int64)
        renum2[util] = np.arange(len(util))
        xy, q = xyz[util].copy(), renum2[q]
        p = xy[q]
        aire = 0.5 * np.sum(p[:, :, 0] * np.roll(p[:, :, 1], -1, axis=1)
                            - np.roll(p[:, :, 0], -1, axis=1) * p[:, :, 1], axis=1)
        q[aire < 0] = q[aire < 0][:, ::-1]

        # lissage laplacien de la zone non structurée (noeuds intérieurs)
        xy = lissage(xy, q, N_LISSAGE, figes=np.unique(q[:n_struct]))
        return xy, q, est_pale

    # La zone non structurée dépend de l'algorithme : on en essaie plusieurs
    # et on garde le maillage dont le pire quadrangle est le meilleur.
    meilleur, qual = None, -np.inf
    for algo in (6, 5, 8, 1):       # frontal-Delaunay, Delaunay, frontal quads, MeshAdapt
        gmsh.model.mesh.clear()
        gmsh.option.setNumber("Mesh.Algorithm", algo)
        try:
            gmsh.model.mesh.generate(2)
            res = extraire()
        except Exception:
            continue
        sj = jacobien_quads(res[0], res[1]).min()
        if sj > qual:
            meilleur, qual = res, sj
        if qual > 0.5:
            break
    gmsh.model.remove()
    if meilleur is None or qual < 0.2:
        raise RuntimeError(f"Empreinte 2D de qualité insuffisante "
                           f"(jacobien normalisé min = {qual:.2f}) : "
                           "ajuster LC_PALE / LC_JANTE.")
    xy, q, est_pale = meilleur
    return xy, q, est_pale, pas, tg


def lissage(xy, q, n_iter, figes=()):
    """Lissage laplacien d'un maillage de quads ; le bord et les noeuds
    'figes' ne bougent pas."""
    xy = xy.copy()
    ar = np.stack([q, np.roll(q, -1, axis=1)], axis=-1).reshape(-1, 2)
    cle, nb_vu = np.unique(np.sort(ar, axis=1), axis=0, return_counts=True)
    fixe = np.zeros(len(xy), dtype=bool)
    fixe[np.unique(cle[nb_vu == 1])] = True                # bord du secteur
    fixe[np.asarray(figes, dtype=np.int64)] = True
    libre = ~fixe
    i, j = np.r_[cle[:, 0], cle[:, 1]], np.r_[cle[:, 1], cle[:, 0]]
    deg = np.bincount(i, minlength=len(xy))[:, None]
    for _ in range(n_iter):
        moy = np.column_stack([np.bincount(i, xy[j, c], len(xy)) for c in (0, 1)])
        xy[libre] = 0.5 * xy[libre] + 0.5 * (moy / deg)[libre]
    return xy


def jacobien_quads(xy, q):
    """Jacobien normalisé (sinus de l'angle) minimal aux 4 coins de chaque quad."""
    p = xy[q]
    e1, e0 = np.roll(p, -1, axis=1) - p, np.roll(p, 1, axis=1) - p
    cr = e1[:, :, 0] * e0[:, :, 1] - e1[:, :, 1] * e0[:, :, 0]
    return (cr / (np.linalg.norm(e1, axis=2) * np.linalg.norm(e0, axis=2))).min(axis=1)


# =============================================================================
# 2. Lois radiales
# =============================================================================
def largeur_disque(r):
    """Largeur axiale du disque w(r) : moyeu / voile / jante, raccords en cosinus."""
    def rampe(r, r0, r1, w0, w1):
        t = np.clip((r - r0) / (r1 - r0), 0.0, 1.0)
        return w0 + (w1 - w0) * 0.5 * (1.0 - np.cos(np.pi * t))
    r = np.asarray(r, dtype=float)
    return np.where(r <= R_VOILE_INT,
                    rampe(r, R_MOYEU, R_VOILE_INT, W_MOYEU, W_VOILE),
                    rampe(r, R_VOILE_EXT, R_JANTE_INT, W_VOILE, W_JANTE))


def rayons_disque():
    """Niveaux radiaux du disque, de l'alésage à la jante : couches fines sous
    le pied de pale (LC_RADIAL_PIED), progression géométrique vers
    LC_RADIAL_DISQUE en descendant."""
    h, d = LC_RADIAL_PIED, [0.0]
    total = R_JANTE - R_ALESAGE
    while d[-1] + 1.5 * h < total:
        d.append(d[-1] + h)
        h = min(h * RAISON_RADIALE, LC_RADIAL_DISQUE)
    d.append(total)
    return R_JANTE - np.array(d[::-1])


# =============================================================================
# 3. Construction du maillage 3D (numpy)
# =============================================================================
def construire():
    xy, q2, est_pale, pas, tg = empreinte_2d()
    n2 = len(xy)
    s2, z2 = xy[:, 0], xy[:, 1]

    # --- appariement des côtés périodiques du secteur ---
    u = s2 - z2 * tg
    tol = 1e-6 * pas
    gauche = np.where(np.abs(u + 0.5 * pas) < tol)[0]
    droite = np.where(np.abs(u - 0.5 * pas) < tol)[0]
    gauche = gauche[np.argsort(z2[gauche])]
    droite = droite[np.argsort(z2[droite])]
    assert len(gauche) == len(droite) and len(gauche) > 1
    assert np.allclose(z2[gauche], z2[droite], atol=1e-6 * W_JANTE), \
        "côtés périodiques non conformes"

    # --- numérotation globale des noeuds du disque : (secteur, noeud 2D, niveau)
    rd = rayons_disque()
    nd = len(rd)                                  # niveaux radiaux du disque
    propre = np.ones(n2, dtype=bool)
    propre[droite] = False                        # appartiennent au secteur suivant
    n_prop = int(propre.sum())
    loc = np.full(n2, -1, dtype=np.int64)
    loc[propre] = np.arange(n_prop)
    N = N_PALES
    ks = np.arange(N)[:, None]
    idx2 = np.empty((N, n2), dtype=np.int64)      # index "colonne" globale
    idx2[:, propre] = ks * n_prop + loc[propre][None, :]
    idx2[:, droite] = ((ks + 1) % N) * n_prop + loc[gauche][None, :]
    n_col = N * n_prop
    # noeud disque = colonne * nd + niveau
    # Sous la jante, l'empreinte est progressivement "détendue" (mêmes
    # éléments, noeuds mieux répartis) : les mailles fines du pied de pale ne
    # sont pas propagées telles quelles jusqu'à l'alésage.
    xy_det = xy
    for n_it in (DETENTE_ITER, DETENTE_ITER // 2, DETENTE_ITER // 4):
        essai = lissage(xy, q2, n_it)
        if jacobien_quads(essai, q2).min() > 0.3:
            xy_det = essai
            break
    b = np.clip((R_JANTE - rd) / (R_JANTE - R_VOILE_EXT), 0.0, 1.0)
    b = b * b * (3.0 - 2.0 * b)                   # 0 à la jante -> 1 au voile
    s_niv = (1.0 - b) * s2[propre][:, None] + b * xy_det[propre, 0][:, None]
    z_niv = (1.0 - b) * z2[propre][:, None] + b * xy_det[propre, 1][:, None]
    ech = largeur_disque(rd) / W_JANTE            # mise à l'échelle axiale
    X = np.empty((N, n_prop, nd, 3))
    for k in range(N):
        th = s_niv / R_JANTE + 2 * np.pi * k / N
        X[k, :, :, 0] = rd * np.cos(th)
        X[k, :, :, 1] = rd * np.sin(th)
        X[k, :, :, 2] = z_niv * ech
    coords = [X.reshape(-1, 3)]
    n_disque = n_col * nd

    def nd_id(k, i, j):                           # i : noeuds 2D (tableau)
        return idx2[k, i] * nd + j

    # --- noeuds de pale (niveaux 1..nb au-dessus de la jante) ---
    nb = N_COUCHES_PALE
    np_pale = np.unique(q2[est_pale])             # noeuds 2D du profil
    assert not np.isin(np_pale, np.concatenate([gauche, droite])).any()
    locp = np.full(n2, -1, dtype=np.int64)
    locp[np_pale] = np.arange(len(np_pale))
    npp = len(np_pale)
    t = np.arange(1, nb + 1) / nb
    phi = np.radians(VRILLAGE_DEG) * t
    k_ech = 1.0 - (1.0 - EFFILEMENT) * t
    rp = R_JANTE + HAUTEUR * t
    sp, zp = s2[np_pale][:, None], z2[np_pale][:, None]
    s_t = k_ech * (sp * np.cos(phi) + zp * np.sin(phi))      # (npp, nb)
    z_t = k_ech * (-sp * np.sin(phi) + zp * np.cos(phi))
    for k in range(N):
        th = s_t / rp + 2 * np.pi * k / N         # épaisseur physique conservée
        Xp = np.stack([rp * np.cos(th), rp * np.sin(th), z_t], axis=-1)
        coords.append(Xp.reshape(-1, 3))
    coords = np.vstack(coords)

    def pale_id(k, i, j):                         # j = 0 : noeud de jante partagé
        if j == 0:
            return nd_id(k, i, nd - 1)
        return n_disque + k * npp * nb + locp[i] * nb + (j - 1)

    # --- hexaèdres ---
    hex_d, hex_p = [], []
    qp = q2[est_pale]
    for k in range(N):
        for j in range(nd - 1):
            hex_d.append(np.hstack([nd_id(k, q2, j), nd_id(k, q2, j + 1)]))
        for j in range(nb):
            hex_p.append(np.hstack([pale_id(k, qp, j), pale_id(k, qp, j + 1)]))
    hex_d, hex_p = np.vstack(hex_d), np.vstack(hex_p)

    # --- peaux (quads, normale sortante) ---
    # arêtes orientées de chaque quad 2D : (a -> b), intérieur à gauche
    ar = np.stack([q2, np.roll(q2, -1, axis=1)], axis=-1).reshape(-1, 2)
    prop = np.repeat(np.arange(len(q2)), 4)
    cle = np.sort(ar, axis=1)
    cle = cle[:, 0] * n2 + cle[:, 1]
    ordre = np.argsort(cle, kind="stable")
    cs = cle[ordre]
    deb = np.r_[True, cs[1:] != cs[:-1]]
    fin = np.r_[cs[1:] != cs[:-1], True]
    bord = ordre[deb & fin]                       # arêtes vues une seule fois
    zm = z2[ar[bord]].mean(axis=1)
    e_amont = ar[bord[zm < -0.5 * W_JANTE + 1e-6]]
    e_aval = ar[bord[zm > +0.5 * W_JANTE - 1e-6]]
    # contour du profil : arêtes partagées entre un quad pale et un quad hors pale
    i1, i2 = ordre[deb & ~fin], ordre[np.r_[False, (deb & ~fin)[:-1]]]
    mixte = est_pale[prop[i1]] != est_pale[prop[i2]]
    e_prof = ar[np.where(est_pale[prop[i1]], i1, i2)[mixte]]

    peau = {g: [] for g in ("ALESAGE", "AMONT", "AVAL", "JANTE",
                            "PEAU_PALES", "TETE_PALES")}
    qh = q2[~est_pale]
    for k in range(N):
        peau["ALESAGE"].append(nd_id(k, q2[:, ::-1], 0))
        peau["JANTE"].append(nd_id(k, qh, nd - 1))
        peau["TETE_PALES"].append(pale_id(k, qp, nb))
        for j in range(nd - 1):
            for nom, e in (("AMONT", e_amont), ("AVAL", e_aval)):
                a, b = e[:, 0], e[:, 1]
                peau[nom].append(np.column_stack(
                    [nd_id(k, a, j), nd_id(k, b, j),
                     nd_id(k, b, j + 1), nd_id(k, a, j + 1)]))
        a, b = e_prof[:, 0], e_prof[:, 1]
        for j in range(nb):
            peau["PEAU_PALES"].append(np.column_stack(
                [pale_id(k, a, j), pale_id(k, b, j),
                 pale_id(k, b, j + 1), pale_id(k, a, j + 1)]))
    peau = {g: np.vstack(v) for g, v in peau.items()}
    return coords, {"DISQUE": hex_d, "PALES": hex_p}, peau


def volumes_hexa(coords, h):
    """Jacobien au centre de chaque hexa (doit être > 0)."""
    p = coords[h]
    d1 = (p[:, [1, 2, 6, 5]] - p[:, [0, 3, 7, 4]]).mean(axis=1)
    d2 = (p[:, [3, 2, 6, 7]] - p[:, [0, 1, 5, 4]]).mean(axis=1)
    d3 = (p[:, 4:] - p[:, :4]).mean(axis=1)
    return np.einsum("ij,ij->i", np.cross(d1, d2), d3)


def elancement(coords, h):
    """Rapport plus grande / plus petite arête de chaque hexa."""
    p = coords[h]
    a = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
         (0, 4), (1, 5), (2, 6), (3, 7)]
    L = np.stack([np.linalg.norm(p[:, i] - p[:, j], axis=1) for i, j in a], 1)
    return L.max(axis=1) / L.min(axis=1)


# =============================================================================
# 4. Injection dans Gmsh et écriture MED
# =============================================================================
def main():
    gui = "--gui" in sys.argv
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    coords, vol, peau = construire()

    for nom, h in vol.items():
        v = volumes_hexa(coords, h)
        assert v.min() > 0, f"hexaèdres retournés dans {nom}"

    gmsh.model.add("roue_aubagee")
    ent = {}
    for dim, groupes in ((3, vol), (2, peau)):
        for nom in groupes:
            ent[nom] = gmsh.model.addDiscreteEntity(dim)
    premier = ent["DISQUE"]
    gmsh.model.mesh.addNodes(3, premier, np.arange(1, len(coords) + 1),
                             coords.ravel())
    for nom, h in vol.items():
        gmsh.model.mesh.addElementsByType(ent[nom], 5, [], (h + 1).ravel())
        gmsh.model.addPhysicalGroup(3, [ent[nom]], name=nom)
    for nom, q in peau.items():
        gmsh.model.mesh.addElementsByType(ent[nom], 3, [], (q + 1).ravel())
        gmsh.model.addPhysicalGroup(2, [ent[nom]], name=nom)

    # noeud maître sur l'axe : maille POI1 seule dans le groupe CENTRE
    n_centre = len(coords) + 1
    ent_c = gmsh.model.addDiscreteEntity(0)
    gmsh.model.mesh.addNodes(0, ent_c, [n_centre], [0.0, 0.0, Z_CENTRE])
    gmsh.model.mesh.addElementsByType(ent_c, 15, [], [n_centre])
    gmsh.model.addPhysicalGroup(0, [ent_c], name="CENTRE")

    n_hex = sum(len(h) for h in vol.values())
    print(f"Noeuds : {len(coords)} + 1 (CENTRE)")
    for nom, h in vol.items():
        el = gmsh.model.mesh.getElementsByType(5, ent[nom])[0]
        sj = gmsh.model.mesh.getElementQualities(el, "minSJ")
        print(f"  {nom:<7}: {len(h):>7} HEXA8   jacobien normalisé min "
              f"{sj.min():.2f} / moy {sj.mean():.2f}   élancement max "
              f"{elancement(coords, h).max():.1f}")
    print(f"Total   : {n_hex} HEXA8, "
          f"{sum(len(q) for q in peau.values())} QUAD4 de peau")

    gmsh.write(FICHIER)
    print(f"Écrit : {FICHIER}")
    if gui:
        gmsh.fltk.run()
    gmsh.finalize()


if __name__ == "__main__":
    main()
