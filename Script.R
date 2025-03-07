# Pacotes necessários
install.packages("AGHmatrix")
install.packages("sommer")
install.packages("rBayesianOptimization")
install.packages("dplyr")
install.packages("ggplot2")

# Observação:
# Este script realiza o cálculo da matriz H no contexto do SSGBLUP. Para executar a análise, é necessário ter previamente calculadas as matrizes de parentesco:
# - Matriz A: Matriz de parentesco baseada em pedigree.
# - Matriz G: Matriz de parentesco genômico, derivada de marcadores moleculares.

# Carregar os pacotes
library(AGHmatrix); library(sommer) library(rBayesianOptimization); library(dplyr); library(ggplot2)

################################ Grid Search ################################

# Definir a grade de parâmetros
(tau_values <- seq(0.1, 2.0, by = 0.1))  # tau de 0.1 a 2.0, com incrementos de 0.1
(omega_values <- seq(-1.0, 1.0, by = 0.1))   # omega de -1.0 a 1.0, com incrementos de 0.1

# Lista para armazenar as matrizes H
H_matrizes <- list()

# Criar uma lista para armazenar os resultados
results <- data.frame()

# Iterar sobre as combinações de tau e omega
for (tau in tau_values) {
  for (omega in omega_values) {
    # Calcular a matriz H com os parâmetros tau e omega
    H <- Hmatrix(A = matriz_A, G = matriz_G, tau = tau, omega = omega, method = "Martini")
    
    # Armazenar a matriz H na lista
    H_matrizes[[paste("tau", tau, "omega", omega, sep = "_")]] <- H        # 420 combinações de tau e omega - ok 
    
    # Ajustar o modelo misto usando a matriz H
    model <- mmer(Y ~ 1, random = ~ vsr(Geno, Gu = H), rcov = ~ units, data = pheno_df)
    
    # Obter os valores preditos
    predicted_values <- predict(model, D = "Geno")$pval
    blups <- predicted_values$predicted.value
    
    # Calcular a capacidade preditiva (correlação entre observado e predito)
    predictive_ability <- cor(pheno_df$Y, blups[match(pheno_df$Geno, predicted_values$Geno)])
    
    # Calcular a inflação (coeficiente de regressão)
    inflation <- coef(lm(pheno_df$Y ~ blups[match(pheno_df$Geno, predicted_values$Geno)]))[2]
    
    # Calcular o Erro Quadrático Médio (MSE) para cada combinação de tau e omega
    mse <- mean((pheno_df$Y - blups[match(pheno_df$Geno, predicted_values$Geno)])^2)
    
    # Armazenar os resultados
    results <- rbind(results, data.frame(
      tau = tau,
      omega = omega,
      predictive_ability = predictive_ability,
      inflation = inflation,
      mse = mse
    ))
  }
}

# Adicionar uma coluna "is_best" para indicar as melhores combinações
# Definir um limite para considerar uma combinação como "melhor"
threshold_predictive_ability <- 0.95 * max(results$predictive_ability)  # 95% da capacidade preditiva máxima
threshold_inflation <- 1.1  # Inflação máxima de 1.1

# Criar a coluna "is_best"
results$is_best <- results$predictive_ability >= threshold_predictive_ability & results$inflation <= threshold_inflation
head(results); tail(results) 

# Filtrar apenas as combinações onde is_best é TRUE
(best_combinations <- subset(results, is_best == TRUE))

# Selecionar a combinação com a maior capacidade preditiva
(best_combination <- best_combinations[which.max(best_combinations$predictive_ability), ])

### Fim da análise 
## Agora que temos a melhor combinação de parâmetros, podemos calcular a matriz H

# Calcular a matriz H com a inclusão de pesos (parametros)
(matriz_H=Hmatrix(matriz_A2, matriz_G, tau=1.9, omega=0.1, method="Martini"))

#-----------------------------------------------------------------------------------------#

################################  Otimização bayesiana ################################ 

# Criar uma lista para armazenar os resultados
results_bayes <- data.frame()

# Função objetivo para otimização bayesiana
objective_function <- function(tau, omega) {
  # Calcular a matriz H com os parâmetros tau e omega
  H <- Hmatrix(A = matriz_A2, G = matriz_G, tau = tau, omega = omega, method = "Martini")
  
  # Ajustar o modelo misto usando a matriz H
  model <- mmer(Y ~ 1, random = ~ vsr(Geno, Gu = H), rcov = ~ units, data = pheno_df)
  
  # Obter os valores preditos
  predicted_values <- predict(model, D = "Geno")$pval
  blups <- predicted_values$predicted.value
  
  # Calcular a capacidade preditiva (correlação entre observado e predito)
  predictive_ability <- cor(pheno_df$Y, blups[match(pheno_df$Geno, predicted_values$Geno)])
  
  # Calcular a inflação (coeficiente de regressão)
  inflation <- coef(lm(pheno_df$Y ~ blups[match(pheno_df$Geno, predicted_values$Geno)]))[2]
  
  # Calcular o Erro Quadrático Médio (MSE)
  mse <- mean((pheno_df$Y - blups[match(pheno_df$Geno, predicted_values$Geno)])^2)
  
  # Armazenar os resultados na tabela
  results_bayes <<- rbind(results_bayes, data.frame(
    tau = tau,
    omega = omega,
    predictive_ability = predictive_ability,
    inflation = inflation,
    mse = mse,  
    row.names = NULL  
  ))
  
  # Retornar a capacidade preditiva como métrica a ser maximizada
  return(list(Score = predictive_ability, Pred = inflation))
}

# Definir os limites para tau e omega
bounds <- list(tau = c(0.1, 2.0),         # tau entre 0.1 e 2.0
               omega = c(-1.0, 1.0))      # omega entre -1.0 e 1.0

# Executar a otimização bayesiana
set.seed(123)                # Definir semente para reproducibilidade
bayes_opt <- BayesianOptimization(
  FUN = objective_function,  # Função objetivo
  bounds = bounds,           # Limites dos parâmetros
  init_points = 10,          # Número de pontos iniciais aleatórios para explorar o espaço de busca
  n_iter = 20,               # Número de iterações após os pontos iniciais
  acq = "ucb",               # Função de aquisição (Upper Confidence Bound)
  kappa = 2.576,             # Parâmetro de exploração (default) - Valores maiores exploram mais
  verbose = TRUE             # Mostrar progresso
)

# Resultados da otimização
bayes_opt

# Extrair o histórico de iterações
(history_bayes <- bayes_opt$History)

# tabela de resultados
(results_bayes_rounded <- results_bayes %>%
  mutate(
    tau = round(tau, 1),
    omega = round(omega, 1))
)

# Extrair os melhores parâmetros
(best_tau <- round(bayes_opt$Best_Par["tau"], 1))        
(best_omega <- round(bayes_opt$Best_Par["omega"], 1))  

# Extrair a capacidade preditiva e a inflação
(best_predictive_ability <- bayes_opt$Best_Value)
(best_inflation <- bayes_opt$Pred)

### Fim da análise 

# Calcular a matriz H com os melhores parâmetros
(matriz_H <- Hmatrix(matriz_A2, matriz_G, tau = best_tau, omega = best_omega, method = "Martini"))

#-----------------------------------------------------------------------------------------#

############# Comparação de resultados das duas abordagens #############

# Melhores resultados do Grid Search
(best_grid_search <- results[which.max(results$predictive_ability), ])   # tau = 1.9; omega = 0.1

# Melhores resultados da Otimização Bayesiana
(best_bayes_opt <- history_bayes[which.max(history_bayes$Value), ])      # tau = 1.8; omega = 0.1

#-----------------------------------------------------------------------------------------#

############## Visualização gráfica dos resultados ##############

#### Grid Search ####

# Gráfico de calor para a capacidade preditiva
ggplot(results, aes(x = tau, y = omega, fill = predictive_ability)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Capacidade Preditiva em Função de Tau e Omega",
       x = "Tau", y = "Omega", fill = "Capacidade Preditiva")

# Gráfico de dispersão 
ggplot(results, aes(x = tau, y = omega, color = predictive_ability)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Grid Search: Capacidade Preditiva",
       x = "Tau", y = "Omega", color = "Capacidade Preditiva")

# Gráfico de inflação
ggplot(results, aes(x = tau, y = omega, fill = inflation)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Inflação em Função de Tau e Omega",
       x = "Tau", y = "Omega", fill = "Inflação")


# Gráfico de calor das melhores combinações 
ggplot(best_combinations, aes(x = tau, y = omega, fill = predictive_ability)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Melhores Combinações de Tau e Omega",
       x = "Tau", y = "Omega", fill = "Capacidade Preditiva")


#### Otimização Bayesiana ####

# Gráfico de calor para a capacidade preditiva
ggplot(results_bayes_rounded, aes(x = tau, y = omega, fill = predictive_ability)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Otimização Bayesiana: Capacidade Preditiva em Função de Tau e Omega",
       x = "Tau", y = "Omega", fill = "Capacidade Preditiva") 

# Gráfico de dispersão 
ggplot(history_bayes, aes(x = tau, y = omega, color = Value)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Otimização Bayesiana: Capacidade Preditiva",
       x = "Tau", y = "Omega", color = "Capacidade Preditiva") 

# Gráfico de calor para a inflação
ggplot(results_bayes_rounded, aes(x = tau, y = omega, fill = inflation)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Otimização Bayesiana: Inflação em Função de Tau e Omega",
       x = "Tau", y = "Omega", fill = "Inflação") 
