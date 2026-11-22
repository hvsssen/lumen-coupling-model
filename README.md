# Lumen coupling model

Modélisation de deux lumens reliés par un flux de soluté, sans pompe osmotique. Le flux est entièrement piloté par la différence de pression de Laplace entre les deux lumens.

## Modèle physique

Chaque lumen est modélisé comme une lentille formée de deux segments circulaires symétriques de part et d'autre de la jonction cellule-cellule, d'où l'aire :

$$A(R) = R^2 \left(2\theta - \sin(2\theta)\right)$$

L'angle $\theta$ ne dépend que des tensions de surface (fixes) :

$$\theta_i = \arccos\left(\frac{\gamma_C}{2\gamma_i}\right)$$

donc il reste constant au cours du temps. En dérivant la relation $A(R)$ et en l'injectant dans le bilan de flux $J_{12} = C_{12}\left(\frac{\gamma_1}{R_1} - \frac{\gamma_2}{R_2}\right)$, on obtient le système d'équations différentielles résolu dans `lumen_two_odes.py`, en $R_1$ et $R_2$ :

$$\frac{dR_1}{dt} = -\frac{J_{12}}{2 g_1 R_1}, \qquad \frac{dR_2}{dt} = \frac{J_{12}}{2 g_2 R_2}$$

où $C_{12}$ est la conductivité de la jonction entre les deux lumens, $\gamma_1$, $\gamma_2$ les tensions de surface des lumens, et $g_1$, $g_2$ les constantes géométriques $g_i = 2\theta_i - \sin(2\theta_i)$ issues de $\theta_1$, $\theta_2$.

Comme `dR1` et `dR2` viennent du même flux `J12` avec des signes opposés, l'aire totale `A1+A2` est conservée : pas de pompe osmotique, donc pas de création ni de destruction de matière, seulement un transfert entre les deux lumens. Le système est instable (type mûrissement d'Ostwald) : le lumen le plus "défavorisé" (rayon plus petit et/ou tension plus forte) se vide dans l'autre.

## Utilisation

```
python3 lumen_two_odes.py
```

Dépendances : `numpy`, `scipy`, `matplotlib`.

Les paramètres (`C12`, `gamma1`, `gamma2`, `gammaC`, `R1_0`, `R2_0`, `t_max`) sont à régler en haut du fichier. Le script affiche l'instant où un lumen se vide (le cas échéant), vérifie la conservation de l'aire totale, et trace R1(t), R2(t), A1(t), A2(t), P1(t), P2(t).

## Suite

Prochaine étape : intégrer la pression osmotique au modèle, en plus de la pression de Laplace déjà prise en compte.
