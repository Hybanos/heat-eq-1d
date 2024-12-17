# Notations:
Pour les questions suivantes on a:
 - scalaires lettres grecques
 - vecteurs lettres latines minuscule
 - matrices lettres latines majuscule

## Exercice 3

---
#### 1:
Comme un pointeur vers double (`double *` et non `double **`)
Par example pour une matrice n*m:
```C
double a[5][3] = {1,1,1,2,3,4,3,5,2,4,2,5,5,4,3};
```
ou:
```C
double *a = malloc(sizeof(double) * n * m)
```

---
#### 2:
C'est une macro qui permet d'indiquer à une fonction si une matrice est sockée par lignes ou par colones.

---
#### 3:
La leading dimension correspond a la taille des segments dans lesquels les éléments de la matrice sont continus en mémoire

---
#### 4:
Multiplication matrice-vecteur sur matrice bande.
$\alpha Ax + \beta y$

---
#### 5:
Factorisation LU de matrice bande en utilisant un pivot.
Replace les valeurs de la matrice A par LU.
Implémentation de BLAS3

---
#### 6:
Résolution d'un système linéaire sur une matrice factorisée LU.
$A\times X = B$

---
#### 7:
Résolution d'un système linéaire. 
Implémente DGBTRF et DGBTRS

---
#### 8:
formule du résidu:
$$r = b - A\hat{x}$$
norme résiduelle arrière:
$$\frac{||b - A\hat{x}||}{||A||||\hat{x}||}$$
norme résiduelle avant:
$$\frac{||r||}{||b||} = \frac{||b - A\hat{x}||}{||b||}$$

pour calculer la norme d'un vecteur `x`:
`n = sqrt(ddot(&n, x, 1, x, 1));`

pour $||x - \hat{x}||$:
`daxpy(&n, -1, x^, 1, x, 1);`


## Exercice 4:

### 2-3:
Soit un vecteur $x$ qu'aquedes1:
on a $AB x = y$ tel que y est un vecteur de nul avec des 1 à chaque extrémités 
$$
\begin{bmatrix}
 0 &-1 &-1 &-1 &-1  \\
 2 & 2 & 2 & 2 & 2  \\
-1 &-1 &-1 &-1 & 0 
\end{bmatrix}
\times
\begin{bmatrix}
1 \\
1 \\
1 \\
1 \\
1 
\end{bmatrix}=
\begin{bmatrix}
1 \\
0 \\
0 \\
0 \\
1 
\end{bmatrix}
$$

