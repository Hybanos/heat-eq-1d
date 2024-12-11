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
$$r = \frac{||x - \hat{x}||}{||x||}$$

## Exercice 4:
