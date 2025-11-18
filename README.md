# Résumé théorique de l’ADPLS

On considère deux matrices décrivant les mêmes \(n\) individus :

- \(X = [x_1,\ldots,x_p]\) : variables numériques centrées  
- \(Y = [y_1,\ldots,y_q]\) : indicatrices (non centrées) d’une variable qualitative à \(q\) modalités  

Les individus sont pondérés par la matrice diagonale \(W = \mathrm{diag}(w_i)\).  
L’espace des variables est muni de la métrique ACP \(M\).

L’Analyse Discriminante PLS (ADPLS) construit des composantes \(f = X M u\) en maximisant le critère :

\[
\|f\|_W^2 \, R^2(f,Y)
\quad\text{sous la contrainte}\quad
u' M u = 1.
\]

On montre que :

\[
\|f\|_W^2\, R^2(f,Y) = \| \hat{X} M u \|_W^2,
\qquad \hat{X} = \Pi_Y X.
\]

---

## 1. Programme de rang 1

Le premier axe discriminant est solution du programme :

\[
(P):\quad \max_{u' M u = 1}\, \|\hat{X} M u\|_W^2,
\quad\text{où}\quad
E = \hat{X}' W \hat{X}.
\]

Ce problème équivaut à :

\[
\max_{u' M u = 1} u' M E M u.
\]

Le lagrangien conduit à l’équation propre :

\[
E M u = \lambda u.
\]

Ainsi, **la première composante** est :

\[
f^1 = X M u_1,
\]

où \(u_1\) est le vecteur propre \(M\)-unitaire associé à la plus grande valeur propre de \(EM\).

---

## 2. Programme de rang \(h\)

Pour les composantes suivantes, on impose l’orthogonalité :

\[
F_{h-1}' W f_h = 0,
\quad
F_{h-1} = [f^1,\ldots,f^{h-1}].
\]

Le programme devient :

\[
(P_2):\quad \max_{\substack{u' M u = 1\\ D' M u = 0}} u' M E M u,
\qquad
D' = F_{h-1}' W X.
\]

Les conditions du premier ordre mènent à :

\[
\Pi_{D^\perp} E M u = \lambda u,
\]

où

\[
\Pi_{D^\perp} = I - D(D' M D)^{-1} D' M
\]

est le projecteur sur l’espace orthogonal à \(\langle D \rangle\).  
La composante s’en déduit : \(f_h = X M u_h\).

---

## 3. Indicateurs d'interprétation

### Inertie expliquée
\[
S(f) = \frac{\|f\|_W^2}{\operatorname{tr}(X' W X)}.
\]

Il s’agit de la part d’inertie totale expliquée par la composante discriminante.

### Pouvoir discriminant
\[
R^2(f,Y) = \frac{\|\hat{X} M u\|_W^2}{\|f\|_W^2}.
\]

- \(R^2 = 1\) : séparation parfaite des classes  
- \(R^2 = 0\) : aucune discrimination

---

## 4. Représentations graphiques

### Coordonnées des individus
Dans le plan \((h,m)\), l’individu \(i\) est représenté par :

\[
(\tilde{f}_{ih},\ \tilde{f}_{im}) ,
\qquad
\tilde{f}_h = \frac{f_h}{\|f_h\|_W}.
\]

### Centres de classes
Avec la matrice indicatrice \(Y\) :

\[
\bar{F}_H = (Y' W Y)^{-1} Y' W F_H.
\]

### Variables dans le plan dual
La variable \(x_j\) est représentée par ses corrélations avec les composantes :

\[
\left(
\frac{\langle x_j , f_h \rangle_W}{\|x_j\|_W \|f_h\|_W},
\;
\frac{\langle x_j , f_m \rangle_W}{\|x_j\|_W \|f_m\|_W}
\right).
\]

---

## 5. Données simulées utilisées dans les exemples

Les données considérées sont une matrice \(X \in \mathbb{R}^{9 \times 7}\),  
une variable qualitative \(Y\in\{A,B,C,D,E\}\),  
et une pondération uniforme \(W = \frac{1}{n} I_n\).

La métrique est :

\[
M = \mathrm{diag}\left(\frac{1}{\sigma_j^2}\right),
\qquad
\sigma_j^2 = x_j' W x_j.
\]

Ce choix suit les recommandations du cours d’Analyse Discriminante :  
prendre \(M = I\) si les variables sont centrées-réduites,  
ou \(M = \mathrm{diag}(\sigma_j^{-2})\) si elles ne sont que centrées.

---

## 6. Code

L’ensemble des fonctions et scripts associés aux calculs et graphiques est présenté dans le dossier \R du package adPLS.
