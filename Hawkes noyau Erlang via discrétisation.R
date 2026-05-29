# Donnees
mu = 1
alpha = 0.5
beta = 1
tau = numeric(0)
T_max = 150
k = 4
delta = 0.1

nb_candidats = 0
nb_acceptes  = 0

temps_pic = (k - 1) / beta

grille_t = seq(0, T_max, by = delta)

# Boucle
for (i in 1:(length(grille_t) - 1)) {
  t_debut = grille_t[i]
  t_fin   = grille_t[i] + delta
  
  if (length(tau) == 0) {
    lambda_bar = mu
  } else {
    dt_min = t_debut - tau
    dt_max = t_fin - tau
    
    noyau_min = (alpha * beta^k * dt_min^(k-1) * exp(-beta * dt_min)) / factorial(k-1)
    noyau_max = (alpha * beta^k * dt_max^(k-1) * exp(-beta * dt_max)) / factorial(k-1)
    
    if (k == 1) {
      sommet = 0
    } else {
      sommet = (alpha * beta * (k-1)^(k-1) * exp(1 - k)) / factorial(k-1)
    }

    max_par_point = ifelse(dt_min <= temps_pic & temps_pic <= dt_max, sommet, pmax(noyau_min, noyau_max))
    
    lambda_bar = mu + sum(max_par_point)
  }
  n_candidats_local = rpois(1, lambda_bar * delta)
  
  if (n_candidats_local > 0) {
    temps_candidats = sort(runif(n_candidats_local, t_debut, t_fin))

    for (t_cand in temps_candidats) {
      nb_candidats = nb_candidats + 1
      dt_candidat = t_cand - tau
      lambda_actual = mu + alpha * sum((beta^k * dt_candidat^(k-1) * exp(-beta * dt_candidat)) / factorial(k-1))
      
      if (runif(1, 0, 1) <= lambda_actual / lambda_bar) {
        tau = c(tau, t_cand)
        nb_acceptes = nb_acceptes + 1
      }
    }
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

