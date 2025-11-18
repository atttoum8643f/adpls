On considère deux matrices de données décrivant $n$ individus (en lignes) :  
$X = [x_1, \ldots, x_p]$ codant $p$ variables numériques centrées, et $Y = [y_1, \ldots, y_q]$ codant une variable qualitative à $q$ modalités, représentée par ses indicatrices non centrées.  

On note $W = \mathrm{diag}(w_i;\, i = 1, \ldots, n)$ la matrice diagonale des poids des individus.  

L’espace $\mathbb{R}^p$ est muni de la métrique d’ACP associée à $X$, notée $M$.  

On rappelle que l’Analyse Discriminante PLS (ADPLS) est fondée sur la maximisation, pour chaque composante $f = X M u$ du critère produit suivant : $\|f\|_W^2 \, R^2(f, Y)$ sous la contrainte $u' M u = 1$.

Montrons que $\|f\|_W^2 \, R^2(f, Y) = \| \hat{X}Mu\|_W^2$, où $\hat{X}= \Pi_Y X$.

Par définition $R^2(f, Y) = \frac{\|\Pi_Y f\|_W^2}{\|f\|_W^2}$ c'est-à-dire:

\begin{align*}
\|f\|_W^2 \, R^2(f, Y) &= \|\Pi_Y f\|_W^2\\
        &= \|\Pi_Y XMu\|_W^2 \\
        &= \|\hat{X}Mu\|_W^2 
\end{align*}

## \textcolor{blue}{1.1. Programme 1 - Rang 1}

Dans ce premier programme, nous cherchons à déterminer la première composante de l’analyse à partir des matrices de données disponibles. 
Il s’agit d’introduire la démarche permettant d’identifier cette première direction principale avant de généraliser au cas de rang supérieur.

Soit $E = \hat{X}'W\hat{X}$ la matrice de variance empirique entre classes. Considérons le programme:

$$
(P): \quad \max_{\;u'Mu = 1} \| \hat{X}Mu \|_W^2
$$

dont la solution fornit la première composante principale $f^1$.

Montrons que $f^1 = XMu_1$, où $u_1$ est le vecteur propre $M$-unitaire de $EM$ associé à la plus grande valeur propre.

On peut réécrire le programme $(P)$:

\begin{align*}
(P): &\quad \max_{\;u'Mu = 1} \| \hat{X}Mu \|_W^2 \\
\iff (P): &\quad \max_{\;u'Mu = 1} u'M\hat{X}'W\hat{X}Mu\\
\iff (P): &\quad \max_{\;u'Mu = 1} u'MEMu
\end{align*}

Le lagrangien associé au problème $(P)$ s'écrit:

$$
\mathcal{L} = u' M E M u - \lambda (u' M u - 1)
$$

D'après la condition du premier ordre on a: $\frac{\partial \mathcal{L}}{\partial \lambda} = 0 \iff u' M u = 1$ et 

\begin{align*}
\frac{\partial \mathcal{L}}{\partial u} = 0 \iff MEMu &= \lambda Mu \quad (1)\\
\iff EMu &= \lambda u \quad (2) \quad \text{(car M est définie positive donc inversible)}
\end{align*}

Ainsi, $u$ est valeur propres de $EM$ associée à la valeur propre $\lambda$, de plus, $u'(1) \iff u' M E M u = \lambda \quad \text{( la valeur propre maximale)} \iff \|\hat{X}Mu_1\|_W^2 = \lambda_1$. Or, $\hat{X}Mu_1 = \Pi_Y XMu_1 = \Pi_Y f_1$, on en déduit alors que $f^1 = XMu_1$, où $u_1$ est le vecteur propre $M$-unitaire de $EM$ associé à la plus grande valeur propre.

Si on pose $u^* = M^{1/2}u, \qquad X^* = X M^{1/2}, \qquad \hat{X}^* = \Pi_Y X^*$, alors:

\begin{align*}
(P): &\quad \max_{\; u'Mu = 1}
\| \hat{X} M u \|_W^2 \\
\iff (P): &\quad \max_{\; u' M^{1/2} M^{1/2}u = 1}
\| \hat{X} M u \|_W^2 \\
\iff (P): &\quad  \max_{\; u^{*'} u^* = 1} 
\| \Pi_Y X M^{1/2} M^{1/2} u \|_W^2 \\
\iff (P): &\quad  \max_{\; u^{*'} u^* = 1} 
\| \hat{X}^* u^* \|_W^2
\end{align*}

$$
\mathcal{L} = u^{*'} E^* u^* - \lambda (u^{*'}u^* - 1)
$$

D'après la condition du premier ordre, on a :  
$\dfrac{\partial \mathcal{L}}{\partial \lambda} = 0 \iff u^{*'} u^* = 1$ et

\begin{align*}
\dfrac{\partial \mathcal{L}}{\partial u^*} = 0 &\iff E^* u^* = \lambda u^* \quad (3)
\end{align*}

Ainsi, $u^*$ est vecteur propre de $E^*$ associé à la valeur propre $\lambda$.  
De plus, en multipliant (3) à gauche par $u^{*'}$, on obtient $u^{*'} E^* u^* = \lambda$, la valeur propre maximale.  

On en déduit alors que $f^1 = XMu_1 = XM^{1/2} M^{1/2}u_1 = X^* u_1^*$.  

Montrons à présent que $u_1^*$ est le vecteur propre $I$-unitaire de $E^*$ associé à la plus grande valeur propre.

Soient $u_i^*, u_j^*$ deux vecteurs propres de la matrice $E^*$ associées aux valeurs propres respectives $\lambda_i$ et $\lambda_j$. Comme $u_i^*$ et $u_j^*$ sont solution du problème:

$$
(Q): \quad  \max_{u^{*'} u^* = 1} 
\| \hat{X}^* u^* \|_W^2
$$

On a alors: $\|  u_i^* \|_I = \|  u_j^* \|_I = 1$. Il suffit de montrer que $\langle u_i^*, u_j^* \rangle = 0$ afin d'obtenir l'orthonormalité des $(u_k^*)_{k=1,\cdots,p}$.

On sait que:

$$
\begin{cases}
E^* u_i^* = \lambda_i u_i^* \quad (4), \\
E^* u_j^* = \lambda_j u_j^* \quad (5).
\end{cases}
$$

c'est-à-dire $u_j^{*'}(4)$ équivaut: 


\begin{align*}
u_j^{*'} E^* u_i^* &= \lambda_i\, u_j^{*'} u_i^* \\
\iff (E^* u_j^*)' u_i^* &= \lambda_i\, u_j^{*'} u_i^* \quad \text{(car $E^*$ est symétrique)} \\
\iff \lambda_j\, u_j^{*'} u_i^* &= \lambda_i\, u_j^{*'} u_i^* \quad \text{(d'après (5))} \\
\iff (\lambda_j - \lambda_i)\, u_j^{*'} u_i^* &= 0
\end{align*}


Or, $\lambda_i \neq \lambda_j$ donc $u_j^{*'} u_i^* = 0$. Ainsi $u_i^*$ et $u_j^*$ sont orthogonaux, et la normalisation précédente donne l'orthonormalité recherchée.

## \textcolor{blue}{1.2. Programme 2 - Rang h}

Dans cette partie, on désire obtenir des composantes deux à deux orthogonales.  
On note $F^{h-1} = [\, f^1, \dots, f^{h-1} \,]$.

La $h$-ième composante $f_h$ doit alors satisfaire la contrainte d’orthogonalité suivante : $F_{h-1}'\, W\, f_h = 0.$

Montrons que le programme :

$$
(P_2): \quad \max_{\substack{u'Mu = 1 \\ D'Mu = 0}} u'MEMu \quad , \quad \text{où} \quad D' = F^{h-1'} W X
$$

nous mène à la recherche de $u$ solution de :

$$
\Pi_D^{\perp} EMu = \lambda u \quad (6), \quad \text{où} \quad \Pi_D^{\perp} = I - D(D'MD)^{-1}D'M \quad \text{et} \quad \lambda \text{ maximale.}
$$

Considérons le lagrangien associé au programme $(P_2)$ :

$$
\mathcal{L} = u'MEMu - \lambda (u'Mu - 1) - \mu' (D'Mu)
$$

Nous avons, par la condition du premier ordre :

$$
\begin{cases}
\frac{\partial \mathcal{L}}{\partial \lambda} = 0 \iff u'Mu = 1, \\
\frac{\partial \mathcal{L}}{\partial \mu} = 0 \iff D'Mu = 0
\end{cases}
$$

Quant à la dérivée de $\mathcal{L}$ selon $u$, nous obtenons :

\begin{align*}
\frac{\partial \mathcal{L}}{\partial u} = 0 &\iff MEMu - \lambda Mu - MD\mu = 0 \quad (i)\\
&\iff EMu - \lambda u - D\mu = 0 \quad (ii)
\end{align*}

En prémultipliant l'équation $(i)$ à gauche par $u'$, nous obtenons :

$$
u'(i) \iff u'MEMu = \lambda \quad \text{(maximale)}
$$

En effet, $u'MD = (D'Mu)' = 0$, car par hypothèse $D'Mu = 0$.

\begin{align*}
D'(i) &\iff D'MEMu = D'MD \mu \\
&\iff  \mu = (D'MD)^{-1} D'MEMu
\end{align*}

On remplace $\mu$ dans $(ii)$, on a :

\begin{align*}
EMu - D(D'MD)^{-1} D'MEMu &= \lambda u\\
\iff (I - D(D'MD)^{-1} D'M) EMu &= \lambda u\\
\iff \Pi_{D^{\perp}} EMu &= \lambda u  \quad (6)
\end{align*}

L'équation $(6)$ peut être réécris de la manière suivante: $\Pi_{D^{\perp}} EM \Pi_{D^{\perp}}u = \lambda u$ à condition que $u$ soit dans l'espace engendré par $D^{\perp}$, c'est-à-dire  $u \in \langle D^{\perp} \rangle$.

Montrons que $u \in \langle D^{\perp} \rangle$.

\begin{align*}
u \in \langle D^{\perp} \rangle &\iff \Pi_{D^{\perp}} u = u\\
&\iff (I - D(D'MD)^{-1} D'M) u =  u\\
&\iff u - D(D'MD)^{-1} (D'M u) =  u\\
&\iff u =  u
\end{align*}

En effet, nous avons dans le programme $(P_2)$ la contrainte $D'M u = 0$ donc $D(D'MD)^{-1} (D'M u) = 0$. Ainsi, $u \in \langle D^{\perp} \rangle$, d'où $\Pi_{D^{\perp}} EM \Pi_{D^{\perp}}u = \lambda u \quad (7)$. En considérant la transformation étoile que l'on a vue précédament l'équation $(7)$ équivaut à:

\begin{align*}
\Pi_{D^{\perp}} E M \Pi_{D^{\perp}} u &= \lambda u \\
\iff \Pi_{D^{\perp}} E \Pi_{D^{\perp}}' M u &= \lambda u \quad (7') \quad \text{(car $\Pi_{D^{\perp}}$ est $M$-symétrique.)}
\end{align*}

En prémultipliant l'équation $(7')$ à gauche par $M^{1/2}$, nous obtenons :

\begin{align*}
M^{1/2} (7') 
&\iff M^{1/2} \Pi_{D^{\perp}} E \Pi_{D^{\perp}}' M^{1/2} M^{1/2} u = \lambda M^{1/2} u \\
&\iff M^{1/2} \Pi_{D^{\perp}} E \Pi_{D^{\perp}}' M^{1/2} u^{*} = \lambda u^{*}  \quad (8) \\
&\iff \tilde{E} u^{*} = \lambda u^{*}
\end{align*}

qui caractérise la diagonalisation de la matrice $\tilde{E}$ symétrique.

## \textcolor{blue}{1.3. Interprétation des indicateurs}

Dans cette partie, on cherche à interpréter les composantes discriminantes à l’aide de deux indicateurs :  

$$
S(f) = \frac{\|f\|_W^2}{\operatorname{tr}(X'WX)}
$$

où $\|f\|_W^2 = f'Wf$ représente l’inertie portée par la composante $f$,  
et $\operatorname{tr}(X'WX)$ l’inertie totale du nuage des points.  

Ainsi, $S(f)$ mesure la proportion de la variance totale expliquée par la composante discriminante $f$.  
C’est l’analogue du pourcentage d’inertie expliquée en ACP, mais dans un cadre discriminant.  

Le pouvoir discriminant de la composante $f$ est donné par :

$$
R^2(f,Y) = \frac{\|\hat{X} M u\|_W^2}{\|f\|_W^2}
$$

avec $\hat{X} = \Pi_Y X$, la projection de $X$ sur l’espace des classes.  

Le numérateur correspond à l’inertie inter-classes, c’est-à-dire la variance expliquée par les classes,  
et le dénominateur à l’inertie totale de la composante $f$.  

Ainsi, $R^2(f,Y)$ représente la proportion de la variance de $f$ due aux différences entre classes.  

Si $R^2(f,Y) = 1$, la séparation entre les classes est parfaite.  
Si $R^2(f,Y) = 0$, il n’y a aucune discrimination.  

C’est l’équivalent du coefficient de détermination $R^2$, appliqué à la structure de classes.


## \textcolor{blue}{1.4. Représentations graphiques}

Dans le plan discriminant $(h,m)$, l'individu $i$ est représenté par ses coordonnées sur les composantes réduites :

$$
\left( \tilde{f}_{ih},\, \tilde{f}_{im} \right)
$$

où

$$
\tilde{f}_h = \frac{f_h}{\|f_h\|_W}, \quad f_h = X M u_h.
$$



Soit $Y$ la matrice indicatrice des classes (dimension $n \times q$), et $W$ la matrice de pondération.  
Le centre de gravité de la classe $k$ sur l'axe $h$ est :

$$
\bar{f}_{kh}
= 
\frac{\displaystyle \sum_{i \in \text{classe } k} w_i f_{ih}}
     {\displaystyle \sum_{i \in \text{classe } k} w_i}.
$$


\paragraph{Forme matricielle.}

$$
\bar{F}_H = (Y'WY)^{-1} Y'W F_H \quad \text{avec} \quad F_H = [f_1, \ldots, f_H].
$$


Pour les composantes réduites $\tilde{F}_H = [\tilde{f}_1, \ldots, \tilde{f}_H]$, on obtient :

$$
\tilde{\bar{F}}_H = (Y'WY)^{-1} Y'W\, \tilde{F}_H.
$$

Cela fournit les coordonnées barycentriques des centres des classes dans l’espace discriminant.


Dans le plan dual, la variable $x_j$ est représentée par ses corrélations avec les 
composantes discriminantes :

$$
\left(
\frac{\langle x_j \mid f_h \rangle_W}{\|x_j\|_W\,\|f_h\|_W},
\;
\frac{\langle x_j \mid f_m \rangle_W}{\|x_j\|_W\,\|f_m\|_W}
\right).
$$


# \textcolor{red}{2. Programmation}

Nous allons à présent mettre en œuvre les concepts précédemment présentés à travers une programmation.  
L’objectif est d’illustrer les différentes étapes de calcul et de vérifier expérimentalement les résultats théoriques obtenus.


Dans la suite, nous considérons les données simulées suivantes, afin de vérifier la bonne fonctionnalité des programmes développés et d’illustrer les différentes étapes de calcul.


On considère le tableau de données simulé :

$$
X =
\begin{bmatrix}
-2.512 & -0.760 & -3.394 & 3.86 & -0.403 & 3.27 & 5.42 \\
0.444 & 7.934 & 1.861 & -2.03 & 4.349 & -1.04 & -3.50 \\
6.566 & -2.110 & -4.890 & 3.39 & 1.140 & -0.31 & 6.28 \\
-2.393 & 7.359 & -4.388 & -1.53 & -0.313 & 6.31 & 4.33 \\
3.927 & -0.183 & 5.342 & 5.96 & 7.778 & 2.20 & 4.41 \\
-4.702 & 2.972 & 7.903 & 5.58 & -4.502 & 1.62 & -1.47 \\
7.283 & 6.772 & -0.422 & -1.94 & -2.633 & -4.79 & -4.26 \\
-1.748 & -0.445 & -3.713 & 4.43 & -0.220 & 6.33 & 3.73 \\
6.016 & 0.101 & 0.398 & -2.00 & -0.166 & -4.50 & 6.46
\end{bmatrix}.
$$

La variable qualitative associée est :
$Y = [A, \ A, \ B, \ C, \ D, \ B, \ C, \ E, \ A]'$.

La matrice des poids des individus est définie uniformément par : $W = \frac{1}{n} I_n$. 

La métrique dans l’espace des variables est donnée par :

$$
M = \operatorname{diag}\!\left(\frac{1}{\sigma_j^2}\right)_{j=1,\cdots,p}
$$

où $M = I_n$, et où chaque $\sigma_j^2 = x_j' W x_j$ représente la variance pondérée de la variable $x_j$, selon que $X$ est respectivement standardisé ou seulement centré.

Ainsi, conformément au cours \textbf{AMULTIVAR\_L2\_2groupes partie 2 — AD (page 46)},  
\begin{quote}
« La régression linéaire, utilisant $M = (X'WX)^{-1}$, est régularisée par la régression PLS1 qui lui substitue $M = I$ si les $x_j$ sont centrées-réduites, ou $M = \operatorname{diag}(\sigma_j^{-2})_j$ si elles ne sont que centrées.  
De même, pour tenir compte de la force structurelle dans $\langle X \rangle$ en AFD, on munit $\mathbb{R}^p$ de la métrique $M = I$ si les $x_j$ sont centrées-réduites (ou $M = \operatorname{diag}(\sigma_j^{-2})_j$ si elles ne sont que centrées). »
\end{quote}


Le code associé à l’ensemble des programmes utilisés pour les calculs et représentations graphiques est présenté en \textbf{annexe}.

