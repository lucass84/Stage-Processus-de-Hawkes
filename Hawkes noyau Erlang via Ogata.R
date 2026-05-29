# Donnees
mu = 1
alpha = 0.5
beta = 1
tau = numeric(0)
t = 0
T_max = 100
k = 2
nb_candidats = 0
nb_acceptes  = 0

temps_pic = (k - 1) / beta

# Boucle
while (t < T_max) {
  t_depart = t 
  
  if (length(tau) == 0) {
    lambda_bar = mu
  } else {
    dt_depart = t_depart - tau
    terme_gauche = (alpha * beta^k * dt_depart^(k-1) * exp(-beta * dt_depart)) / factorial(k-1)

    if (k == 1) {
      terme_droite = 0
    } else {
      terme_droite = (alpha * beta * (k-1)^(k-1) * exp(1 - k)) / factorial(k-1)
    }
    max_par_point = ifelse(dt_depart > temps_pic, terme_gauche, terme_droite)
    lambda_bar = mu + sum(max_par_point)
  }
  
  U = runif(1, 0, 1)
  w = -log(U) / lambda_bar
  t = t + w 
  
  if (t >= T_max) break
  
  nb_candidats = nb_candidats + 1
  
  if (length(tau) == 0) {
    lambda_actual = mu
  } else {
    dt_candidat = t - tau
    lambda_actual = mu + alpha * sum((beta^k * dt_candidat^(k-1) * exp(-beta * dt_candidat)) / factorial(k-1))
  }
  
  if (runif(1, 0, 1) <= lambda_actual / lambda_bar) {
    tau = c(tau, t)
    nb_acceptes = nb_acceptes + 1
  }
}

prop_acceptes = (nb_acceptes / nb_candidats) * 100
lambda_hat = length(tau) / T_max 

# Affichage
print(lambda_hat)
print(prop_acceptes)
grille_t = seq(0, T_max, length.out = 500)
v_lambda = numeric(length(grille_t))

for (i in 1:length(grille_t)) {
  t_inst = grille_t[i]
  tau_passes = tau[tau < t_inst]
  
  if (length(tau_passes) == 0) {
    v_lambda[i] = mu
  } else {
    v_lambda[i] = mu + alpha * sum((beta^k * (t_inst - tau_passes)^(k-1) * exp(-beta * (t_inst - tau_passes))) / factorial(k-1))
  }
}

plot(grille_t, v_lambda, type = "l", 
     lwd = 2, col = "blue",
     main = "Évolution de l'intensité lambda(t) associée aux sauts",
     xlab = "Temps (t)", 
     ylab = "Intensité lambda(t)", 
     xlim = c(0, T_max))

abline(h = mu, col = "grey", lty = 2)

fonction_escalier = stepfun(tau, 0:length(tau))

plot(fonction_escalier, 
     verticals = TRUE,    
     do.points = FALSE,    
     lwd = 2,
     main = "Processus de comptage cumulé N(t)",
     xlab = "Temps (t)", 
     ylab = "Nombre total d'événements N(t)",
     xlim = c(0, T_max))

# Fonction K de Ripley
lambda_hat = length(tau) / T_max
grille_r = seq(0, T_max / 3, 0.2)
K_hat = numeric(length(grille_r))

for (i in 1:length(grille_r)) {
  r = grille_r[i]
  indices_x = which(tau <= (T_max - r))
  N_x = length(indices_x)
  sum = 0 
  for (j in indices_x) {
    for (k in 1:length(tau)) {
      diff_temps = tau[k] - tau[j]
      if (diff_temps > 0 && diff_temps <= r) {
        sum = sum + 1
      }
    }
    K_hat[i] = (1 / lambda_hat) * (sum / N_x)
  }
}

plot(grille_r, K_hat, type = 'l', main = "Fonction K de Ripley empirique", xlab = "Distance (r)", ylab = "K(r)")