# ADPLS : Analyse Discriminante PLS

Le package **ADPLS** implémente une méthode factorielle discriminante inspirée des approches PLS.  
Il permet d’extraire des **axes discriminants**, de représenter les **individus**, **centres de classes**, et de produire un **cercle des corrélations** des variables.

---

# Sommaire

- [Présentation](#présentation)  
- [Principe de la méthode](#principe-de-la-méthode)  
- [Composantes discriminantes](#composantes-discriminantes)  
- [Indicateurs](#indicateurs)  
- [Graphiques](#graphiques)  
- [Exemple d’utilisation](#exemple-dutilisation)  
- [Annexe mathématique](#annexe-mathématique)

---

# Présentation

Étant donnés :

- une matrice de variables numériques  
  \( X = [x_1,\dots,x_p] \), de dimension \( n \times p \)
- une variable qualitative à \( q \) modalités, codée par  
  \( Y = [y_1,\dots,y_q] \)
- une matrice de pondération diagonale  
  \( W = \mathrm{diag}(w_1,\dots,w_n) \)

L’objectif d’ADPLS est d’extraire des **composantes discriminantes**  
\[
f_h = X M u_h,
\]
qui maximisent la séparation entre groupes.

La matrice \(M\) est une métrique dans l’espace des variables (souvent \(M = I\) pour des données standardisées).

---

# Principe de la méthode

Pour une direction \(u\), on maximise :

\[
\|f\|_W^2\, R^2(f,Y) \quad \text{où } f = XMu.
\]

On montre que :

\[
\|f\|_W^2\,R^2(f,Y)
=
\|\hat{X} M u\|_W^2,
\qquad \hat{X} = \Pi_Y X.
\]

La matrice inter-classes pondérée est :

\[
E = \hat{X}' W \hat{X}.
\]

ADPLS cherche donc des directions maximisant l’inertie inter-classes.

---

# Composantes discriminantes

## Première composante

Le premier axe discriminant est obtenu en résolvant :

\[
\max_{u'Mu = 1} u'MEMu.
\]

Cela conduit à l’équation :

\[
EMu = \lambda u.
\]

La composante associée est :

\[
f_1 = X M u_1.
\]

## Composantes suivantes

On impose l’orthogonalité \(W\)-pondérée :

\[
F_{h-1}'W f_h = 0.
\]

On résout alors :

\[
\max_{\substack{u’Mu = 1 \\ D’Mu = 0}} u'MEMu,
\qquad
D' = F^{h-1'}WX.
\]

Ce qui mène à :

\[
\Pi_D^{\perp} EMu = \lambda u.
\]

---

# Indicateurs

## Inertie expliquée

\[
S(f)
=
\frac{\|f\|_W^2}{\mathrm{tr}(X'WX)}.
\]

## Pouvoir discriminant

\[
R^2(f,Y)
=
\frac{\|\hat{X}Mu\|_W^2}{\|f\|_W^2}.
\]

---
