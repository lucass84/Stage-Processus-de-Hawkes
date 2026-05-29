mu =50
alpha = 0.99999
beta = 1
tau = numeric(0)
t = 0
T_max = 50
r_fixe = 1.0
delta = 0.1
points_proposes = 0
points_acceptes = 0

while (t < T_max) {
  lambda_bar = mu + sum(alpha * exp(-beta * (tau[length(tau)] - tau)))
  
  U = runif(1, 0, 1)
  w = -log(U) / lambda_bar
  t = t + w
  
  if (t >= T_max) break
  
  points_proposes = points_proposes + 1
  
  lambda_actual = mu + sum(alpha * exp(-beta * (t - tau)))
  
  if (runif(1, 0, 1) <= lambda_actual / lambda_bar) {
    tau = c(tau, t)
    points_acceptes = points_acceptes + 1
  }
}

lambda_hat = length(tau) / T_max 
taux_acceptation = (points_acceptes / points_proposes) * 100

print(tau)
cat("Intensité moyenne globale (lambda_hat) :", lambda_hat, "\n")
cat("Taux d'acceptation de l'algorithme :", round(taux_acceptation, 2), "%\n")

grille_t = seq(0, T_max, length.out = 500)
v_lambda = numeric(length(grille_t))

for (i in 1:length(grille_t)) {
  t_inst = grille_t[i]
  tau_passes = tau[tau < t_inst]
  if (length(tau_passes) == 0) {
    v_lambda[i] = mu
  } else {
    v_lambda[i] = mu + alpha * sum(exp(-beta * (t_inst - tau_passes)))
  }
}

plot(grille_t, v_lambda, type = "l", 
     lwd = 2, col = "blue",
     main = "Évolution de l'intensité lambda(t)",
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