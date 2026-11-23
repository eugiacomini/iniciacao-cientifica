library(glmnet) #Lasso e Enet
library(mvnfast) #Normal multivariada menos lenta
library(HDInterval) #Intervalos de Credibilidade HDI
library(compiler) #Compilador que acelera loops
library(profvis) #Reconhecedor de gargalos
library(ggplot2) #Gráficos em geral
library(corrplot) #Gráfico de Correlação
library(MASS) #Modelos StepWise/Backward/Foward
#library(qgraph) #Gráficos em rede
library(openxlsx)
library(readxl)

#Gerando os dados do modelo Y = Xβ + erro
gera_dados = function(tamanho_amostra, sigmaqua_erro, vetor_parametros) {
  n = tamanho_amostra
  p = length(vetor_parametros)
  
  #Declarando o vetor de erros como normal, a verossimilhança como normal e o vetor de betas
  erro = rnorm(n, 0, sqrt(sigmaqua_erro))
  
  matriz_x = matrix(0, nrow = n, ncol = p)
  matriz_x[,1] = 1
  
  for(i in 2:p){
    matriz_x[,i] = rnorm(n, 0, 1)
  }
  
  #Construindo o vetor y resposta
  vetor_y = matriz_x %*% vetor_parametros + erro
  
  return(data.frame("y" = vetor_y, matriz_x))
}

#Modelo de Regressão por Mínimos Quadrados
minimos_quadrados = function(vetor_y, matriz_x) {
  if(nrow(matriz_x)<=ncol(matriz_x)){return(list(rep(0, times = ncol(matriz_x)), 0, rep(0, times = ncol(matriz_x)), "Modelo impossível p>=n"))}
  
  x = matriz_x
  p = length(x[1,])
  n = length(x[,1])
  
  modelo = lm(vetor_y~., data = as.data.frame(x[,-1]))
  sigmaquaErro = (1 / (n - p)) * t(vetor_y) %*% (diag(1, nrow = n, ncol = n) - x %*% solve(t(x) %*% x) %*% t(x)) %*% vetor_y
  beta = as.vector(modelo$coefficients)
  indicadora = c(1, ((anova(modelo)$`Pr(>F)`<0.05)*1)[-ncol(x)])
  
  return(list(beta, sigmaquaErro, indicadora, summary(modelo)))
}

#Modelo de seleção por Single Marker Analysis
sma = function(vetor_y, matriz_x){
  alf_sig = 0.05
  
  p_valores = vector()
  p_valores[1] = 0
  indicadora = vector()
  indicadora[1] = 1
  for(i in 2:ncol(matriz_x)){
    modelo = lm(vetor_y~matriz_x[,i])
    p_valores[i] = summary(modelo)$coefficients[,4][2]
  }
  
  #Correção de bonferroni, selecionando variáveis com p-valor menor que alf_sig/p
  indicadora = ifelse(p_valores > alf_sig/ncol(matriz_x), 0, 1)
  indicadora[1] = 1
  
  beta = matrix(0, nrow = ncol(matriz_x), ncol = 1)
  
  if(sum(indicadora==1)<=length(vetor_y)){
    beta[indicadora==1] = lm(vetor_y~matriz_x[,which(indicadora==1)]-1)$coefficients
  }
  
  residuos = vetor_y - (matriz_x %*% beta)
  sigmaqua_erro = sum(residuos^2) / (length(vetor_y) - sum(indicadora))
  
  return(list(beta, sigmaqua_erro, indicadora))
}

#Modelo de seleção por Backward
backward = function(vetor_y, matriz_x){
  if(nrow(matriz_x)<=ncol(matriz_x)){return(list(rep(0, times = ncol(matriz_x)), 0, rep(0, times = ncol(matriz_x)), "Modelo impossível p>=n"))}
  
  df_x = as.data.frame(matriz_x)
  modelo_inicial = lm(vetor_y~.-1, data = df_x)
  modelo_backward = stepAIC(modelo_inicial, direction = "backward", trace = FALSE)
  indicadora = as.integer(colnames(df_x) %in% names(coef(modelo_backward)))
  indicadora[1] = 1
  
  beta = rep(0, ncol(matriz_x))
  nomes_variaveis = names(coef(modelo_backward))
  
  for (i in 1:length(nomes_variaveis)) {
    indice = which(colnames(df_x) == nomes_variaveis[i])
    if (length(indice) > 0) {
      beta[indice] = coef(modelo_backward)[nomes_variaveis[i]]
    }
  }
  
  sigmaqua_erro = sum((vetor_y - matriz_x %*% beta)^2) / (length(vetor_y) - sum(indicadora))
  
  return(list(beta, sigmaqua_erro, indicadora, summary(modelo_backward)))
}

#Modelo de seleção por Forward
forward = function(vetor_y, matriz_x){
  if(nrow(matriz_x)<=ncol(matriz_x)){return(list(rep(0, times = ncol(matriz_x)), 0, rep(0, times = ncol(matriz_x)), "Modelo impossível p>=n"))}
  
  df_x = as.data.frame(matriz_x)
  names(df_x)[1] = "x0"
  
  modelo_inicial = lm(vetor_y~.-1, data = df_x)
  modelo_nulo = lm(y_treino~x0-1, data = df_x)
  modelo_forward = stepAIC(modelo_nulo, scope = list(lower = modelo_nulo, upper = modelo_inicial), direction = "forward", trace = FALSE)
  indicadora = as.integer(colnames(df_x) %in% names(coef(modelo_forward)))
  indicadora[1] = 1
  
  beta = rep(0, ncol(matriz_x))  
  nomes_variaveis = names(coef(modelo_forward))
  
  for (i in 1:length(nomes_variaveis)) {
    indice = which(colnames(df_x) == nomes_variaveis[i])
    if (length(indice) > 0) {
      beta[indice] = coef(modelo_forward)[nomes_variaveis[i]]
    }
  }
  
  sigmaqua_erro = sum((vetor_y - matriz_x %*% beta)^2) / (length(vetor_y) - sum(indicadora))
  
  return(list(beta, sigmaqua_erro, indicadora, summary(modelo_forward)))
}

#Modelo de seleção por StepWise
stepwise = function(vetor_y, matriz_x){
  if(nrow(matriz_x)<=ncol(matriz_x)){return(list(rep(0, times = ncol(matriz_x)), 0, rep(0, times = ncol(matriz_x)), "Modelo impossível p>=n"))}
  
  df_x = as.data.frame(matriz_x)
  modelo_inicial = lm(vetor_y~.-1, data = df_x)
  modelo_stepwise = stepAIC(modelo_inicial, direction = "both", trace = F)
  indicadora = as.integer(colnames(df_x) %in% names(coef(modelo_stepwise)))
  indicadora[1] = 1
  
  beta = rep(0, ncol(matriz_x))  
  nomes_variaveis = names(coef(modelo_stepwise))
  if (is.null(nomes_variaveis)) {
    nomes_variaveis = character(0)
  }
  
  for (i in 1:length(nomes_variaveis)) {
    indice = which(colnames(df_x) == nomes_variaveis[i])
    if (length(indice) > 0) {
      beta[indice] = coef(modelo_stepwise)[nomes_variaveis[i]]
    }
  }
  
  sigmaqua_erro = sum((vetor_y - matriz_x %*% beta)^2) / (length(vetor_y) - sum(indicadora))
  
  return(list(beta, sigmaqua_erro, indicadora, summary(modelo_stepwise)))
}

#Modelo de seleção por LASSO
lasso = function(vetor_y, matriz_x){
  #Ajuste do modelo LASSO com validação cruzada
  cv_lasso = cv.glmnet(matriz_x[,-1], vetor_y, alpha = 1)
  melhor_lambda = cv_lasso$lambda.min
  lm_lasso = glmnet(matriz_x[,-1], vetor_y, alpha = 1, lambda = melhor_lambda)
  
  #Extração dos coeficientes beta (vetor numérico)
  beta = as.vector(coef(lm_lasso, s = melhor_lambda))
  
  #Indicadora de variáveis selecionadas
  indicadora = ifelse(beta == 0, 0, 1)
  indicadora[1] = 1
  
  #Cálculo da variância do erro (sigma^2)
  sigmaqua_erro = sum((vetor_y - (matriz_x %*% beta))^2) / (length(vetor_y) - sum(indicadora))
  
  return(list(beta, sigmaqua_erro, indicadora))
}

#Modelo de seleção por Elastic Net
enet = function(vetor_y, matriz_x){
  #Ajuste do modelo Elastic Net com validação cruzada
  cv_enet = cv.glmnet(matriz_x[,-1], vetor_y, alpha = 0.5)
  melhor_lambda = cv_enet$lambda.min
  lm_lasso = glmnet(matriz_x[,-1], vetor_y, alpha = 0.5, lambda = melhor_lambda)
  
  #Extração dos coeficientes beta (vetor numérico)
  beta = as.vector(coef(lm_lasso, s = melhor_lambda))
  
  #Indicadora de variáveis selecionadas
  indicadora = ifelse(beta == 0, 0, 1)
  indicadora[1] = 1
  
  #Cálculo da variância do erro (sigma^2)
  sigmaqua_erro = sum((vetor_y - (matriz_x %*% beta))^2) / (length(vetor_y) - sum(indicadora))
  
  return(list(beta, sigmaqua_erro, indicadora))
}

#Modelo de Regressão por MCMC - Gibbs Sampler
gibbs_sampler = function(vetor_y, matriz_x, iter, sigmaqua_beta) {
  
  #Declarando as variáveis pra deixar resumido
  p = ncol(matriz_x)
  a = 0.01
  b = 0.01
  y = vetor_y
  x = matriz_x
  n_amostra = length(y)
  
  #Declarando o vetor da variância do erro e a matriz, que serão preenchidos pelas estimativas na iteração
  sigmaqua_erro = numeric(iter)
  beta = matrix(0, nrow = p, ncol = iter)
  
  #Dando o chute inicial para o amostrador de Gibbs, com um vetor nulo mi_beta de dimensão p
  mi_beta = numeric(p)
  
  #Núcleo do amostrador de Gibbs, iterando sobre uma normal p-variada que irá gerar vetores de betas e uma gama-inversa que gerará valores de sigmaqua_erro
  enableJIT(3)
  for(i in 2:iter) {
    sigmaqua_erro[i] = 1 / rgamma(1, a + n_amostra / 2, b + (1 / 2) * t(y - (x %*% beta[, i - 1])) %*% (y - (x %*% beta[, i - 1])))
    
    matriz_teste = matrix(0, nrow = p, ncol = p)
    matriz_teste[cbind(1:p, 1:p)] = 1/sigmaqua_beta
    matriz_V = solve(matriz_teste + t(x) %*% x / sigmaqua_erro[i])
    matriz_V = (matriz_V + t(matriz_V))/2
    matriz_VA = matriz_V %*% (mi_beta / sigmaqua_beta + t(x) %*% y / sigmaqua_erro[i])
    
    beta[, i] = rmvn(1, mu = matriz_VA, sigma = matriz_V)
  }
  
  #Criando os intervalos de credibilidade para cada beta e selecionando (1) ou excluindo (0) a variável do modelo
  intervalos = matrix(data = 1, nrow = p, ncol = 3)
  enableJIT(3)
  beta = beta[,seq(from=1000,to=iter,by=5)]
  for(j in 1:p){
    intervalos[j,1] = hdi(beta[j,])[1]
    intervalos[j,2] = hdi(beta[j,])[2]
    if(intervalos[j,1]<0 && intervalos[j,2]>0){intervalos[j,3] = 0}
    else{intervalos[j,3] = 1}
  }
  
  return(list(beta, sigmaqua_erro[seq(from=1000,to=iter,by=5)], intervalos))
}

#Modelo de seleção por Spike and Slab
spike_slab = function(vetor_y, matriz_x, iter, sigmaqua_beta) {
  
  #Variáveis suporte e chutes iniciais
  x = matriz_x
  p = length(x[1, ])
  a = 0.01
  b = 0.01
  pro_ber = 0.5
  y = vetor_y
  n_amostra = length(y)
  
  mi_beta = matrix(0, ncol = 1, nrow = p)
  sigmaqua_erro = numeric(iter)
  beta = matrix(0, ncol = iter, nrow = p)
  sigmaqua_erro[1] = 100
  
  indicadora = matrix(1, ncol = iter, nrow = p)
  matriz_sigma = matrix(0, nrow = p, ncol = p)
  matriz_sigma[cbind(1:p, 1:p)] = sigmaqua_beta
  x_reduzida = x
  mi_beta_reduzida = mi_beta
  
  enableJIT(3)
  for (i in 2:iter) {
    ind = indicadora[,i-1]==1
    
    #Declarando matrizes e vetores iniciais
    x_reduzida = x[, ind, drop = FALSE]
    mi_beta_reduzida = mi_beta[ind, drop = FALSE]
    
    mat = t(x_reduzida)%*%x_reduzida/sigmaqua_erro[i-1]
    diag(mat) = diag(mat)+1/sigmaqua_beta
    matriz_V = solve(mat)
    matriz_V = (matriz_V + t(matriz_V))/2
    
    #matriz_V = solve(diag(1/diag(matriz_sigma[ind, ind, drop = FALSE])) + t(x_reduzida) %*% x_reduzida / sigmaqua_erro[i - 1])
    matriz_A = mi_beta_reduzida / sigmaqua_beta + t(x_reduzida) %*% y / sigmaqua_erro[i - 1]
    
    #Gerando os betas
    beta[ind, i] = rmvn(1, matriz_V%*%matriz_A, matriz_V)
    beta[!ind, i] = 0
    
    #Gerando a variância do erro
    sigmaqua_erro[i] = 1 / rgamma(1, a + n_amostra / 2, b + (1 / 2) * t(y - (x %*% beta[, i])) %*% (y - (x %*% beta[, i])))
    
    indicadora[, i] = indicadora[, i - 1]
    
    for (j in 2:p) {
      #Calculando a probabilidade incluindo a variável
      gamma = indicadora[, i]
      gamma[j] = 1
      ind_gamma = gamma==1
      p_gamma = sum(ind_gamma)
      
      x_reduzida = x[, ind_gamma, drop = FALSE]
      mi_beta_reduzida = mi_beta[ind_gamma, drop = FALSE]
      mat = t(x_reduzida)%*%x_reduzida/sigmaqua_erro[i]
      diag(mat) = diag(mat)+1/sigmaqua_beta
      matriz_V = solve(mat)
      matriz_V = (matriz_V + t(matriz_V))/2
      matriz_A = mi_beta_reduzida / sigmaqua_beta + t(x_reduzida) %*% y / sigmaqua_erro[i]
      
      post_ber1 = log(pro_ber) + (determinant(matriz_V, logarithm = TRUE)$modulus)/2 - p_gamma*log(sigmaqua_beta)/2 - (sum(mi_beta_reduzida^2)*sigmaqua_beta - t(matriz_A)%*%matriz_V%*%matriz_A)/2
      
      #Calculando a probabilidade sem incluir a variável
      gamma = indicadora[, i]
      gamma[j] = 0
      ind_gamma = gamma==1
      p_gamma = sum(ind_gamma)
      
      x_reduzida = x[, ind_gamma, drop = FALSE]
      mi_beta_reduzida = mi_beta[ind_gamma, drop = FALSE]
      mat = t(x_reduzida)%*%x_reduzida/sigmaqua_erro[i]
      diag(mat) = diag(mat)+1/sigmaqua_beta
      matriz_V = solve(mat)
      matriz_V = (matriz_V + t(matriz_V))/2
      matriz_A = mi_beta_reduzida / sigmaqua_beta + t(x_reduzida) %*% y / sigmaqua_erro[i]
      
      post_ber0 = log(1 - pro_ber) + (determinant(matriz_V, logarithm = TRUE)$modulus)/2 - p_gamma*log(sigmaqua_beta)/2 - (sum(mi_beta_reduzida^2)*sigmaqua_beta - t(matriz_A)%*%matriz_V%*%matriz_A)/2
      
      #Calculando a probabilidade própria e gerando um valor aleatório da bernoulli com a probabilidade calculada
      maximo = max(post_ber1, post_ber0)
      indicadora[j, i] = rbinom(1, 1, exp(post_ber1 - maximo) / (exp(post_ber1 - maximo) + exp(post_ber0 - maximo)))
    }
  }
  
  return(list(beta[,seq(from=1000,to=iter,by=5)], sigmaqua_erro[seq(from=1000,to=iter,by=5)], indicadora[,seq(from=1000,to=iter,by=5)]))
}

#Modelo de seleção por Data Driven Reversible Jump
ddrj = function(vetor_y, matriz_x, iter, sigmaqua_beta){
  
  #Variáveis de suporte e chutes iniciais
  a = 0.01
  b = 0.01
  sigmaqua_erro = numeric(iter)
  sigmaqua_erro[1] = 100
  ind = 0
  ind_nova = 0
  ind_nova_nova = 0
  
  y = vetor_y
  x = matriz_x
  n_amostra = length(y)
  p_total = length(matriz_x[1,])
  
  sorteia_variavel = NULL
  indicadora = matrix(1, ncol = iter, nrow = p_total)
  p_modelo = length(indicadora[,1]==1)
  mi_beta = matrix(0, ncol = 1, nrow = p_total)
  beta = matrix(0, ncol = iter, nrow = p_total)
  matriz_sigma = matrix(0, ncol = p_total, nrow = p_total)
  matriz_sigma[cbind(1:p_total, 1:p_total)] = sigmaqua_beta
  
  dnorm_atual = sum(dnorm(y,x%*%beta[,1], sqrt(sigmaqua_erro[1]), log = TRUE)) # mudei aqui o sqrt no sigmaqua_erro
  priori_atual = sum(dnorm(beta[,1], mi_beta, sqrt(sigmaqua_beta), log = TRUE)) + dgamma(1/sigmaqua_erro[1], a, b,log=TRUE) # mudei aqui o sqrt no sigmaqua_beta
  
  enableJIT(3)
  for(i in 2:iter){
    
    #Decisão 1: matar (1) uma covariável ou nascer (0) uma covariável
    if(p_modelo==1){p_nasce_ou_morre = 0}
    if(p_total>p_modelo && p_modelo>1){p_nasce_ou_morre = 0.5}
    if(p_modelo==p_total){p_nasce_ou_morre = 1}
    nasce_ou_morre = rbinom(1, 1, p_nasce_ou_morre)
    
    #Decisão 2: sortear uma covariável para nascer ou morrer
    prob_ida = 0
    
    if(nasce_ou_morre==1){
      filtra = which(indicadora[,i-1]==1)[-1]
      if (length(filtra)==1) sorteia_variavel = filtra else sorteia_variavel = sample(filtra, size = 1) 
      indicadora[,i] = indicadora[,i-1]
      indicadora[sorteia_variavel,i] = 0
      prob_ida = log(p_nasce_ou_morre) -log(length(which(indicadora[,i-1]==1)))
    }
    if(nasce_ou_morre==0){
      filtra = which(indicadora[,i-1]==0)
      if (length(filtra)==1) sorteia_variavel = filtra else sorteia_variavel = sample(filtra, size = 1)
      indicadora[,i] = indicadora[,i-1]
      indicadora[sorteia_variavel,i] = 1
      prob_ida = log(p_nasce_ou_morre) -log(length(which(indicadora[,i-1]==0)))
    }
    
    ind = indicadora[,i]==1
    
    #Decisão 3: Gerando os betas e os sigmas do modelo candidato (novo)
    matriz_sigma_atual = matriz_sigma[ind, ind,drop=FALSE]
    x_atual = x[,ind, drop = FALSE]
    mi_beta_atual = mi_beta[ind, drop = FALSE]
    matriz_V = solve(solve(matriz_sigma_atual) + t(x_atual) %*% x_atual / sigmaqua_erro[i - 1])
    matriz_V = (matriz_V+t(matriz_V))/2
    matriz_A = mi_beta_atual / sigmaqua_beta + t(x_atual) %*% y / sigmaqua_erro[i - 1]
    matriz_VA = matriz_V %*% matriz_A
    
    beta[ind,i] = mvnfast::rmvn(1, mu = matriz_VA, sigma = matriz_V)
    beta[!ind,i] = 0
    media_teste_teste = beta[, i][which(beta[, i] != 0)]
    sigmaqua_erro[i] = 1 / rgamma(1, a + n_amostra / 2, b + (1 / 2) * t(y - (x %*% beta[, i])) %*% (y - (x %*% beta[, i])))
    prob_ida = prob_ida + mvnfast::dmvn(media_teste_teste, mu = matriz_VA, sigma = matriz_V, log = TRUE) + dgamma(1/sigmaqua_erro[i], a + n_amostra / 2, b + (1 / 2) * t(y - (x %*% beta[, i])) %*% (y - (x %*% beta[, i])),log=TRUE)
    
    #Decisão 4: verificar se a razão de verossimilhança entre o modelo novo e o modelo antigo é menor que 1. Se for menor que 1, optamos pelo modelo antigo. Se for maior que 1, mantemos a indicadora do jeito que está.
    dnorm_proposto = sum(dnorm(y, x_atual%*%media_teste_teste, sqrt(sigmaqua_erro[i]), log = TRUE)) # mudei aqui o sqrt
    priori_proposta = sum(dnorm(media_teste_teste, mi_beta_atual, sqrt(sigmaqua_beta), log = TRUE)) + dgamma(1/sigmaqua_erro[i], a, b,log=TRUE) # mudei aqui o sqrt
    
    prob_volta = 0
    p_modelo_proposto=sum(indicadora[,i]==1)
    
    if(p_modelo_proposto==1){
      p_nasce_ou_morre_proposto = 0}
    if(p_total>p_modelo_proposto && p_modelo_proposto>1){
      p_nasce_ou_morre_proposto = 0.5
      prob_volta = prob_volta + log(1-p_nasce_ou_morre_proposto)}
    if(p_modelo_proposto==p_total){
      p_nasce_ou_morre_proposto = 1}
    
    ind_nova = indicadora[,i-1]==1
    
    matriz_sigma_atual = matriz_sigma[ind_nova, ind_nova,drop=FALSE]
    x_atual = x[,ind_nova, drop = FALSE]
    mi_beta_atual = mi_beta[ind_nova, drop = FALSE]
    matriz_V = solve(solve(matriz_sigma_atual) + t(x_atual) %*% x_atual / sigmaqua_erro[i])
    matriz_V = (matriz_V+t(matriz_V))/2
    matriz_A = mi_beta_atual / sigmaqua_beta + t(x_atual) %*% y / sigmaqua_erro[i]
    matriz_VA = matriz_V %*% matriz_A
    
    prob_volta = prob_volta - log(length(which(indicadora[,i]==(1-nasce_ou_morre)))) + mvnfast::dmvn(beta[,i-1][indicadora[,i-1] != 0], mu = matriz_VA, sigma = matriz_V, log = TRUE) + dgamma(1/sigmaqua_erro[i-1], a + n_amostra / 2, b + (1 / 2) * t(y - (x %*% beta[,i-1])) %*% (y - (x %*% beta[,i-1])),log=TRUE) # troquei o beta dentro do beta por indicadora no dmvn
    
    if(rbinom(1, 1, min(1, exp(dnorm_proposto - dnorm_atual + priori_proposta - priori_atual + prob_volta - prob_ida)))==0){
      indicadora[,i]=indicadora[,i-1]
      sigmaqua_erro[i]=sigmaqua_erro[i - 1]}
    
    ind_nova_nova = indicadora[,i]==1
    
    p_modelo = sum(indicadora[,i]==1)
    matriz_sigma_atual = matriz_sigma[ind_nova_nova, ind_nova_nova,drop=FALSE]
    x_atual = x[,ind_nova_nova, drop = FALSE]
    mi_beta_atual = mi_beta[ind_nova_nova, drop = FALSE]
    matriz_V = solve(solve(matriz_sigma_atual) + t(x_atual) %*% x_atual / sigmaqua_erro[i]) #troquei o índice do sigmaqua_erro por causa da mudança que fiz acima
    matriz_V = (matriz_V+t(matriz_V))/2
    matriz_A = mi_beta_atual / sigmaqua_beta + t(x_atual) %*% y / sigmaqua_erro[i] #troquei o índice do sigmaqua_erro por causa da mudança que fiz acima
    matriz_VA = matriz_V%*%matriz_A
    
    beta[ind_nova_nova,i] = mvnfast::rmvn(1, matriz_VA, matriz_V)
    beta[!ind_nova_nova,i] = 0
    
    sigmaqua_erro[i] = 1 / rgamma(1, a + n_amostra / 2, b + (1 / 2) * t(y - (x%*%beta[,i])) %*% (y - (x%*%beta[,i])))
    
    media_teste = beta[, i][which(beta[, i] != 0)]
    dnorm_atual = sum(dnorm(y, x_atual%*%media_teste, sqrt(sigmaqua_erro[i]), log = TRUE)) #troquei aqui o sqrt
    priori_atual = sum(dnorm(media_teste, mi_beta_atual, sqrt(sigmaqua_beta), log = TRUE)) + dgamma(1/sigmaqua_erro[i], a, b,log=TRUE) #troquei aqui o sqrt
  }
  return(list(beta[,seq(from=1000,to=iter,by=5)], sigmaqua_erro[seq(from=1000,to=iter,by=5)], indicadora[,seq(from=1000,to=iter,by=5)]))
}

#Função que retorna as métricas para medir a qualidade dos três métodos
qualidade = function(n_dados, var_erro, parametros, sigmaqua_beta, iteracoes, replicas){
  p = length(parametros)
  gibbs_tempo = vector()
  sslab_tempo = vector()
  ddrj_tempo = vector()
  gibbs_reqmbeta = vector()
  sslab_reqmbeta = vector()
  ddrj_reqmbeta = vector()
  gibbs_sens = vector()
  sslab_sens = vector()
  ddrj_sens = vector()
  gibbs_esp = vector()
  sslab_esp = vector()
  ddrj_esp = vector()
  gibbs_reqmpred = vector()
  sslab_reqmpred = vector()
  ddrj_reqmpred = vector()
  gibbs_absevar = vector()
  sslab_absevar = vector()
  ddrj_absevar = vector()
  
  nome_arquivo = paste0(
    "varBeta_", sigmaqua_beta,
    "_varErro_", var_erro,
    "_nDados_", n_dados,
    ".xlsx"
  )
  
  arquivo = createWorkbook()
  addWorksheet(arquivo, "gibbs_beta")
  addWorksheet(arquivo, "gibbs_sigmaqua")
  addWorksheet(arquivo, "gibbs_intervalos")
  addWorksheet(arquivo, "sslab_beta")
  addWorksheet(arquivo, "sslab_sigmaqua")
  addWorksheet(arquivo, "sslab_indicadora")
  addWorksheet(arquivo, "ddrj_beta")
  addWorksheet(arquivo, "ddrj_sigmaqua")
  addWorksheet(arquivo, "ddrj_indicadora")
  
  for(i in 1:replicas){
    dados = gera_dados(n_dados,var_erro, parametros)
    treino_ind = sample(seq_len(nrow(dados)), size = floor(0.70*nrow(dados)))
    
    x_treino = as.matrix(dados[,-1][treino_ind,])
    y_treino = c(dados[,1][treino_ind])
    
    x_teste = as.matrix(dados[,-1][-treino_ind,])
    y_teste = c(dados[,1][-treino_ind])
    
    gibbs_tempo[i] = system.time(gibbs <- gibbs_sampler(y_treino, x_treino, iteracoes, sigmaqua_beta))[3]
    sslab_tempo[i] = system.time(sslab <- spike_slab(y_treino, x_treino, iteracoes, sigmaqua_beta))[3]
    ddrj_tempo[i] = system.time(ddrj <- ddrj(y_treino, x_treino, iteracoes, sigmaqua_beta))[3]
    
    gibbs_sens[i] = length(intersect(which(parametros!=0), which(gibbs[[3]][,3]==1)))/length(which(parametros!=0))
    sslab_sens[i] = length(intersect(which(parametros!=0), which(rowMeans(sslab[[3]])>0.5)))/length(which(parametros!=0))
    ddrj_sens[i] = length(intersect(which(parametros!=0), which(rowMeans(ddrj[[3]])>0.5)))/length(which(parametros!=0))
    
    gibbs_esp[i] = length(intersect(which(parametros==0), which(gibbs[[3]][,3]==0)))/length(which(parametros==0))
    sslab_esp[i] = length(intersect(which(parametros==0), which(rowMeans(sslab[[3]])<0.5)))/length(which(parametros==0))
    ddrj_esp[i] = length(intersect(which(parametros==0), which(rowMeans(ddrj[[3]])<0.5)))/length(which(parametros==0))
    
    gibbs_reqmbeta[i] = sqrt(sum((rowMeans(gibbs[[1]])-parametros)^2)/length(parametros))
    sslab_reqmbeta[i] = sqrt(sum((rowMeans(sslab[[1]])-parametros)^2)/length(parametros))
    ddrj_reqmbeta[i] = sqrt(sum((rowMeans(ddrj[[1]])-parametros)^2)/length(parametros))
    
    sel_gibbs = which(gibbs[[3]][,3]==1)
    if(length(sel_gibbs)==1){gibbs_reqmpred[i] = sqrt(sum((x_teste[,sel_gibbs]*mean(gibbs[[1]][sel_gibbs,])-y_teste)^2)/length(y_teste))}
    if(length(sel_gibbs)!=1){gibbs_reqmpred[i] = sqrt(sum((x_teste[,sel_gibbs]%*%t(t(rowMeans(gibbs[[1]][sel_gibbs,])))-y_teste)^2)/length(y_teste))}
    
    sel_sslab = which(rowMeans(sslab[[3]])>0.5)
    if(length(sel_sslab)==1){sslab_reqmpred[i] = sqrt(sum((x_teste[,sel_sslab]*mean(sslab[[1]][sel_sslab,])-y_teste)^2)/length(y_teste))}
    if(length(sel_sslab)!=1){sslab_reqmpred[i] = sqrt(sum((x_teste[,sel_sslab]%*%t(t(rowMeans(sslab[[1]][sel_sslab,])))-y_teste)^2)/length(y_teste))}
    
    sel_ddrj = which(rowMeans(ddrj[[3]])>0.5)
    if(length(sel_ddrj)==1){ddrj_reqmpred[i] = sqrt(sum((x_teste[,sel_ddrj]*mean(ddrj[[1]][sel_ddrj,])-y_teste)^2)/length(y_teste))}
    if(length(sel_ddrj)!=1){ddrj_reqmpred[i] = sqrt(sum((x_teste[,sel_ddrj]%*%t(t(rowMeans(ddrj[[1]][sel_ddrj,])))-y_teste)^2)/length(y_teste))}
    
    gibbs_absevar[i] = abs(mean(gibbs[[2]][-c(1, 2, 3)])-var_erro)
    sslab_absevar[i] = abs(mean(sslab[[2]][-c(1, 2, 3)])-var_erro)
    ddrj_absevar[i] = abs(mean(ddrj[[2]][-c(1, 2, 3)])-var_erro)
    
    writeData(arquivo, sheet = "gibbs_beta", x = paste0("Amostra"," ", i), startRow = (i-1)*(p+1)+i)
    writeData(arquivo, sheet = "gibbs_beta", x = gibbs[[1]], startRow = (i-1)*(p+1)+i+1, colNames = FALSE)
    
    writeData(arquivo, sheet = "gibbs_sigmaqua", x = paste0("Amostra"," ", i), startRow = 3*(i-1)+1)
    writeData(arquivo, sheet = "gibbs_sigmaqua", t(gibbs[[2]]), startRow = 3*(i-1)+2, colNames = FALSE)
    
    writeData(arquivo, sheet = "gibbs_intervalos", x = paste0("Amostra"," ", i), startRow = (i-1)*(p+1)+i)
    writeData(arquivo, sheet = "gibbs_intervalos", gibbs[[3]], startRow = (i-1)*(p+1)+i+1, colNames = FALSE)
    
    writeData(arquivo, sheet = "sslab_beta", x = paste0("Amostra"," ", i), startRow = (i-1)*(p+1)+i)
    writeData(arquivo, sheet = "sslab_beta", sslab[[1]], startRow = (i-1)*(p+1)+i+1, colNames = FALSE)
    
    writeData(arquivo, sheet = "sslab_sigmaqua", x = paste0("Amostra"," ", i), startRow = 3*(i-1)+1)
    writeData(arquivo, sheet = "sslab_sigmaqua", t(sslab[[2]]), startRow = 3*(i-1)+2, colNames = FALSE)
    
    writeData(arquivo, sheet = "sslab_indicadora", x = paste0("Amostra"," ", i), startRow = (i-1)*(p+1)+i)
    writeData(arquivo, sheet = "sslab_indicadora", sslab[[3]], startRow = (i-1)*(p+1)+i+1, colNames = FALSE)
    
    writeData(arquivo, sheet = "ddrj_beta", x = paste0("Amostra"," ", i), startRow = (i-1)*(p+1)+i)
    writeData(arquivo, sheet = "ddrj_beta", ddrj[[1]], startRow = (i-1)*(p+1)+i+1, colNames = FALSE)
    
    writeData(arquivo, sheet = "ddrj_sigmaqua", x = paste0("Amostra"," ", i), startRow = 3*(i-1)+1)
    writeData(arquivo, sheet = "ddrj_sigmaqua", t(ddrj[[2]]), startRow = 3*(i-1)+2, colNames = FALSE)
    
    writeData(arquivo, sheet = "ddrj_indicadora", x = paste0("Amostra"," ", i), startRow = (i-1)*(p+1)+i)
    writeData(arquivo, sheet = "ddrj_indicadora", ddrj[[3]], startRow = (i-1)*(p+1)+i+1, colNames = FALSE)
  }
  
  tabela = data.frame("TEMPO_GIBBS" = gibbs_tempo,
                      "TEMPO_SSLAB" = sslab_tempo,
                      "TEMPO_DDRJ" = ddrj_tempo,
                      
                      "SENS_GIBBS" = gibbs_sens,
                      "SENS_SSLAB" = sslab_sens,
                      "SENS_DDRJ" = ddrj_sens,
                      
                      "ESP_GIBBS" = gibbs_esp,
                      "ESP_SSLAB" = sslab_esp,
                      "ESP_DDRJ" = ddrj_esp,
                      
                      "REQMPRED_GIBBS" = gibbs_reqmpred,
                      "REQMPRED_SSLAB" = sslab_reqmpred,
                      "REQMPRED_DDRJ" = ddrj_reqmpred,
                      
                      "REQMBETA_GIBBS" = gibbs_reqmbeta,
                      "REQMBETA_SSLAB" = sslab_reqmbeta,
                      "REQMBETA_DDRJ" = ddrj_reqmbeta,
                      
                      "ABSEVAR_GIBBS" = gibbs_absevar,
                      "ABSEVAR_SSLAB" = sslab_absevar,
                      "ABSEVAR_DDRJ" = ddrj_absevar)
  
  addWorksheet(arquivo, "metricas")
  writeData(arquivo, sheet = "metricas", tabela, rowNames = FALSE, colNames = TRUE)
  saveWorkbook(arquivo, file.path("/home/giacomini/Testes", nome_arquivo), overwrite = TRUE)
  metricas = read_excel(file.path("/home/giacomini/Testes", nome_arquivo), sheet = "metricas")
  #return(View(metricas))
}

#Função que retorna as métricas para medir a qualidade do SMA, LASSO e Elastic Net
qualidade_2 = function(n_dados, var_erro, parametros, replicas){
  p = length(parametros)
  sma_tempo = vector()
  lasso_tempo = vector()
  enet_tempo = vector()
  sma_reqmbeta = vector()
  lasso_reqmbeta = vector()
  enet_reqmbeta = vector()
  sma_sens = vector()
  lasso_sens = vector()
  enet_sens = vector()
  sma_esp = vector()
  lasso_esp = vector()
  enet_esp = vector()
  sma_reqmpred = vector()
  lasso_reqmpred = vector()
  enet_reqmpred = vector()
  sma_absevar = vector()
  lasso_absevar = vector()
  enet_absevar = vector()
  
  nome_arquivo = paste0(
    "varErro_", var_erro,
    "_nDados_", n_dados,
    "_nCovar_", length(parametros),
    ".xlsx"
  )
  
  arquivo = createWorkbook()
  addWorksheet(arquivo, "sma_beta")
  addWorksheet(arquivo, "sma_sigmaqua")
  addWorksheet(arquivo, "sma_intervalos")
  addWorksheet(arquivo, "lasso_beta")
  addWorksheet(arquivo, "lasso_sigmaqua")
  addWorksheet(arquivo, "lasso_indicadora")
  addWorksheet(arquivo, "enet_beta")
  addWorksheet(arquivo, "enet_sigmaqua")
  addWorksheet(arquivo, "enet_indicadora")
  
  for(i in 1:replicas){
    dados = gera_dados(n_dados,var_erro, parametros)
    treino_ind = sample(seq_len(nrow(dados)), size = floor(0.70*nrow(dados)))
    
    x_treino = as.matrix(dados[,-1][treino_ind,])
    y_treino = c(dados[,1][treino_ind])
    
    x_teste = as.matrix(dados[,-1][-treino_ind,])
    y_teste = c(dados[,1][-treino_ind])
    
    sma_tempo[i] = system.time(sma <- sma(y_treino, x_treino))[3]
    lasso_tempo[i] = system.time(lasso <- lasso(y_treino, x_treino))[3]
    enet_tempo[i] = system.time(enet <- enet(y_treino, x_treino))[3]
    
    sma_sens[i] = length(intersect(which(parametros!=0),which(sma[[3]]==1)))/length(which(parametros!=0))
    lasso_sens[i] = length(intersect(which(parametros!=0),which(lasso[[3]]==1)))/length(which(parametros!=0))
    enet_sens[i] = length(intersect(which(parametros!=0), which(enet[[3]]==1)))/length(which(parametros!=0))
    
    sma_esp[i] = length(intersect(which(parametros==0), which(sma[[3]]==0)))/length(which(parametros==0))
    lasso_esp[i] = length(intersect(which(parametros==0), which(lasso[[3]]==0)))/length(which(parametros==0))
    enet_esp[i] = length(intersect(which(parametros==0), which(enet[[3]]==0)))/length(which(parametros==0))
    
    sma_reqmbeta[i] = sqrt(sum((sma[[1]]-parametros)^2)/length(parametros))
    lasso_reqmbeta[i] = sqrt(sum((lasso[[1]]-parametros)^2)/length(parametros))
    enet_reqmbeta[i] = sqrt(sum((enet[[1]]-parametros)^2)/length(parametros))
    
    #Cálculo do REQM para SMA
    sel_sma = which(sma[[3]] == 1)
    if (length(sel_sma) == 1) {
      sma_reqmpred[i] = sqrt(sum((x_teste[, sel_sma] * mean(sma[[1]][sel_sma]) - y_teste)^2) / length(y_teste))
    } else {
      sma_reqmpred[i] = sqrt(sum((x_teste[, sel_sma] %*% sma[[1]][sel_sma] - y_teste)^2) / length(y_teste))
    }
    
    #Cálculo do REQM para LASSO
    sel_lasso = which(lasso[[3]] == 1)
    if (length(sel_lasso) == 1) {
      lasso_reqmpred[i] = sqrt(sum((x_teste[, sel_lasso] * lasso[[1]][sel_lasso] - y_teste)^2) / length(y_teste))
    } else {
      lasso_reqmpred[i] = sqrt(sum((x_teste[, sel_lasso] %*% lasso[[1]][sel_lasso] - y_teste)^2) / length(y_teste))
    }
    
    #Cálculo do REQM para Elastic Net
    sel_enet = which(enet[[3]] == 1)
    if (length(sel_enet) == 1) {
      enet_reqmpred[i] = sqrt(sum((x_teste[, sel_enet] * enet[[1]][sel_enet] - y_teste)^2) / length(y_teste))
    } else {
      enet_reqmpred[i] = sqrt(sum((x_teste[, sel_enet] %*% enet[[1]][sel_enet] - y_teste)^2) / length(y_teste))
    }
    
    sma_absevar[i] = abs(sma[[2]]-var_erro)
    lasso_absevar[i] = abs(lasso[[2]]-var_erro)
    enet_absevar[i] = abs(enet[[2]]-var_erro)
    
    writeData(arquivo, sheet = "sma_beta", x = paste0("Amostra"," ", i), startRow = (i-1)*(p+1)+i)
    writeData(arquivo, sheet = "sma_beta", x = sma[[1]], startRow = (i-1)*(p+1)+i+1, colNames = FALSE)
    
    writeData(arquivo, sheet = "sma_sigmaqua", x = paste0("Amostra"," ", i), startRow = 3*(i-1)+1)
    writeData(arquivo, sheet = "sma_sigmaqua", t(sma[[2]]), startRow = 3*(i-1)+2, colNames = FALSE)
    
    writeData(arquivo, sheet = "sma_intervalos", x = paste0("Amostra"," ", i), startRow = (i-1)*(p+1)+i)
    writeData(arquivo, sheet = "sma_intervalos", sma[[3]], startRow = (i-1)*(p+1)+i+1, colNames = FALSE)
    
    writeData(arquivo, sheet = "lasso_beta", x = paste0("Amostra"," ", i), startRow = (i-1)*(p+1)+i)
    writeData(arquivo, sheet = "lasso_beta", lasso[[1]], startRow = (i-1)*(p+1)+i+1, colNames = FALSE)
    
    writeData(arquivo, sheet = "lasso_sigmaqua", x = paste0("Amostra"," ", i), startRow = 3*(i-1)+1)
    writeData(arquivo, sheet = "lasso_sigmaqua", t(lasso[[2]]), startRow = 3*(i-1)+2, colNames = FALSE)
    
    writeData(arquivo, sheet = "lasso_indicadora", x = paste0("Amostra"," ", i), startRow = (i-1)*(p+1)+i)
    writeData(arquivo, sheet = "lasso_indicadora", lasso[[3]], startRow = (i-1)*(p+1)+i+1, colNames = FALSE)
    
    writeData(arquivo, sheet = "enet_beta", x = paste0("Amostra"," ", i), startRow = (i-1)*(p+1)+i)
    writeData(arquivo, sheet = "enet_beta", enet[[1]], startRow = (i-1)*(p+1)+i+1, colNames = FALSE)
    
    writeData(arquivo, sheet = "enet_sigmaqua", x = paste0("Amostra"," ", i), startRow = 3*(i-1)+1)
    writeData(arquivo, sheet = "enet_sigmaqua", t(enet[[2]]), startRow = 3*(i-1)+2, colNames = FALSE)
    
    writeData(arquivo, sheet = "enet_indicadora", x = paste0("Amostra"," ", i), startRow = (i-1)*(p+1)+i)
    writeData(arquivo, sheet = "enet_indicadora", enet[[3]], startRow = (i-1)*(p+1)+i+1, colNames = FALSE)
  }
  
  tabela = data.frame("TEMPO_sma" = sma_tempo,
                      "TEMPO_lasso" = lasso_tempo,
                      "TEMPO_enet" = enet_tempo,
                      
                      "SENS_sma" = sma_sens,
                      "SENS_lasso" = lasso_sens,
                      "SENS_enet" = enet_sens,
                      
                      "ESP_sma" = sma_esp,
                      "ESP_lasso" = lasso_esp,
                      "ESP_enet" = enet_esp,
                      
                      "REQMPRED_sma" = sma_reqmpred,
                      "REQMPRED_lasso" = lasso_reqmpred,
                      "REQMPRED_enet" = enet_reqmpred,
                      
                      "REQMBETA_sma" = sma_reqmbeta,
                      "REQMBETA_lasso" = lasso_reqmbeta,
                      "REQMBETA_enet" = enet_reqmbeta,
                      
                      "ABSEVAR_sma" = sma_absevar,
                      "ABSEVAR_lasso" = lasso_absevar,
                      "ABSEVAR_enet" = enet_absevar)
  
  addWorksheet(arquivo, "metricas")
  writeData(arquivo, sheet = "metricas", tabela, rowNames = FALSE, colNames = TRUE)
  saveWorkbook(arquivo, file.path("/home/giacomini/Testes", nome_arquivo), overwrite = TRUE)
  metricas = read_excel(file.path("/home/giacomini/Testes", nome_arquivo), sheet = "metricas")
  #return(View(metricas))
}

#Função que retorna as métricas para medir a qualidade do StepWise, Backward e Forward
qualidade_3 = function(n_dados, var_erro, parametros, replicas){
  p = length(parametros)
  stepwise_tempo = vector()
  forward_tempo = vector()
  backward_tempo = vector()
  stepwise_reqmbeta = vector()
  forward_reqmbeta = vector()
  backward_reqmbeta = vector()
  stepwise_sens = vector()
  forward_sens = vector()
  backward_sens = vector()
  stepwise_esp = vector()
  forward_esp = vector()
  backward_esp = vector()
  stepwise_reqmpred = vector()
  forward_reqmpred = vector()
  backward_reqmpred = vector()
  stepwise_absevar = vector()
  forward_absevar = vector()
  backward_absevar = vector()
  
  nome_arquivo = paste0(
    "varErro_", var_erro,
    "_nDados_", n_dados,
    "_nCovar_", length(parametros),
    ".xlsx"
  )
  
  arquivo = createWorkbook()
  addWorksheet(arquivo, "stepwise_beta")
  addWorksheet(arquivo, "stepwise_sigmaqua")
  addWorksheet(arquivo, "stepwise_indicadora")
  addWorksheet(arquivo, "forward_beta")
  addWorksheet(arquivo, "forward_sigmaqua")
  addWorksheet(arquivo, "forward_indicadora")
  addWorksheet(arquivo, "backward_beta")
  addWorksheet(arquivo, "backward_sigmaqua")
  addWorksheet(arquivo, "backward_indicadora")
  
  for(i in 1:replicas){
    dados = gera_dados(n_dados,var_erro, parametros)
    treino_ind = sample(seq_len(nrow(dados)), size = floor(0.70*nrow(dados)))
    
    x_treino = as.matrix(dados[,-1][treino_ind,])
    y_treino = c(dados[,1][treino_ind])
    
    x_teste = as.matrix(dados[,-1][-treino_ind,])
    y_teste = c(dados[,1][-treino_ind])
    
    stepwise_tempo[i] = system.time(stepwise <- stepwise(y_treino, x_treino))[3]
    forward_tempo[i] = system.time(forward <- forward(y_treino, x_treino))[3]
    backward_tempo[i] = system.time(backward <- backward(y_treino, x_treino))[3]
    
    stepwise_sens[i] = length(intersect(which(parametros!=0),which(stepwise[[3]]==1)))/length(which(parametros!=0))
    forward_sens[i] = length(intersect(which(parametros!=0),which(forward[[3]]==1)))/length(which(parametros!=0))
    backward_sens[i] = length(intersect(which(parametros!=0), which(backward[[3]]==1)))/length(which(parametros!=0))
    
    stepwise_esp[i] = length(intersect(which(parametros==0), which(stepwise[[3]]==0)))/length(which(parametros==0))
    forward_esp[i] = length(intersect(which(parametros==0), which(forward[[3]]==0)))/length(which(parametros==0))
    backward_esp[i] = length(intersect(which(parametros==0), which(backward[[3]]==0)))/length(which(parametros==0))
    
    stepwise_reqmbeta[i] = sqrt(sum((stepwise[[1]]-parametros)^2)/length(parametros))
    forward_reqmbeta[i] = sqrt(sum((forward[[1]]-parametros)^2)/length(parametros))
    backward_reqmbeta[i] = sqrt(sum((backward[[1]]-parametros)^2)/length(parametros))
    
    #Cálculo do REQM para stepwise
    sel_stepwise = which(stepwise[[3]] == 1)
    if (length(sel_stepwise) == 1) {
      stepwise_reqmpred[i] = sqrt(sum((x_teste[, sel_stepwise] * mean(stepwise[[1]][sel_stepwise]) - y_teste)^2) / length(y_teste))
    } else {
      stepwise_reqmpred[i] = sqrt(sum((x_teste[, sel_stepwise] %*% stepwise[[1]][sel_stepwise] - y_teste)^2) / length(y_teste))
    }
    
    #Cálculo do REQM para forward
    sel_forward = which(forward[[3]] == 1)
    if (length(sel_forward) == 1) {
      forward_reqmpred[i] = sqrt(sum((x_teste[, sel_forward] * forward[[1]][sel_forward] - y_teste)^2) / length(y_teste))
    } else {
      forward_reqmpred[i] = sqrt(sum((x_teste[, sel_forward] %*% forward[[1]][sel_forward] - y_teste)^2) / length(y_teste))
    }
    
    #Cálculo do REQM para Elastic Net
    sel_backward = which(backward[[3]] == 1)
    if (length(sel_backward) == 1) {
      backward_reqmpred[i] = sqrt(sum((x_teste[, sel_backward] * backward[[1]][sel_backward] - y_teste)^2) / length(y_teste))
    } else {
      backward_reqmpred[i] = sqrt(sum((x_teste[, sel_backward] %*% backward[[1]][sel_backward] - y_teste)^2) / length(y_teste))
    }
    
    stepwise_absevar[i] = abs(stepwise[[2]]-var_erro)
    forward_absevar[i] = abs(forward[[2]]-var_erro)
    backward_absevar[i] = abs(backward[[2]]-var_erro)
    
    writeData(arquivo, sheet = "stepwise_beta", x = paste0("Amostra"," ", i), startCol = i, startRow = 1)
    writeData(arquivo, sheet = "stepwise_beta", x = stepwise[[1]], startCol = i, startRow = 2, colNames = FALSE)
    
    writeData(arquivo, sheet = "stepwise_sigmaqua", x = paste0("Amostra"," ", i), startCol = i, startRow = 1)
    writeData(arquivo, sheet = "stepwise_sigmaqua", t(stepwise[[2]]), startCol = i, startRow = 2, colNames = FALSE)
    
    writeData(arquivo, sheet = "stepwise_indicadora", x = paste0("Amostra"," ", i), startCol = i, startRow = 1)
    writeData(arquivo, sheet = "stepwise_indicadora", stepwise[[3]], startCol = i, startRow = 2, colNames = FALSE)
    
    writeData(arquivo, sheet = "forward_beta", x = paste0("Amostra"," ", i), startCol = i, startRow = 1)
    writeData(arquivo, sheet = "forward_beta", forward[[1]], startCol = i, startRow = 2, colNames = FALSE)
    
    writeData(arquivo, sheet = "forward_sigmaqua", x = paste0("Amostra"," ", i), startCol = i, startRow = 1)
    writeData(arquivo, sheet = "forward_sigmaqua", t(forward[[2]]), startCol = i, startRow = 2, colNames = FALSE)
    
    writeData(arquivo, sheet = "forward_indicadora", x = paste0("Amostra"," ", i), startCol = i, startRow = 1)
    writeData(arquivo, sheet = "forward_indicadora", forward[[3]], startCol = i, startRow = 2, colNames = FALSE)
    
    writeData(arquivo, sheet = "backward_beta", x = paste0("Amostra"," ", i), startCol = i, startRow = 1)
    writeData(arquivo, sheet = "backward_beta", backward[[1]], startCol = i, startRow = 2, colNames = FALSE)
    
    writeData(arquivo, sheet = "backward_sigmaqua", x = paste0("Amostra"," ", i), startCol = i, startRow = 1)
    writeData(arquivo, sheet = "backward_sigmaqua", t(backward[[2]]), startCol = i, startRow = 2, colNames = FALSE)
    
    writeData(arquivo, sheet = "backward_indicadora", x = paste0("Amostra"," ", i), startCol = i, startRow = 1)
    writeData(arquivo, sheet = "backward_indicadora", backward[[3]], startCol = i, startRow = 2, colNames = FALSE)
  }
  
  tabela = data.frame("TEMPO_stepwise" = stepwise_tempo,
                      "TEMPO_forward" = forward_tempo,
                      "TEMPO_backward" = backward_tempo,
                      
                      "SENS_stepwise" = stepwise_sens,
                      "SENS_forward" = forward_sens,
                      "SENS_backward" = backward_sens,
                      
                      "ESP_stepwise" = stepwise_esp,
                      "ESP_forward" = forward_esp,
                      "ESP_backward" = backward_esp,
                      
                      "REQMPRED_stepwise" = stepwise_reqmpred,
                      "REQMPRED_forward" = forward_reqmpred,
                      "REQMPRED_backward" = backward_reqmpred,
                      
                      "REQMBETA_stepwise" = stepwise_reqmbeta,
                      "REQMBETA_forward" = forward_reqmbeta,
                      "REQMBETA_backward" = backward_reqmbeta,
                      
                      "ABSEVAR_stepwise" = stepwise_absevar,
                      "ABSEVAR_forward" = forward_absevar,
                      "ABSEVAR_backward" = backward_absevar)
  
  addWorksheet(arquivo, "metricas")
  writeData(arquivo, sheet = "metricas", tabela, rowNames = FALSE, colNames = TRUE)
  saveWorkbook(arquivo, file.path("/home/giacomini/Testes", nome_arquivo), overwrite = TRUE)
  metricas = read_excel(file.path("/home/giacomini/Testes", nome_arquivo), sheet = "metricas")
  #return(View(metricas))
}

#Análise Descritiva e tratamento de dados
dados = train
nomes_antigos = colnames(train)
names(dados)[names(dados) == "critical_temp"] = "y"
colunas_x = setdiff(colnames(dados), "y")
nomes_novos = paste0("x", seq_along(colunas_x))
names(dados)[names(dados) != "y"] = nomes_novos
dados = dados[, c("y", setdiff(names(dados), "y"))]

#Análise Descritiva
#View(dados)
#View(head(dados))
sum(is.na(dados))
str(dados)
boxplot(scale(dados[,-1]), xaxt = "n", xlab = "Covariáveis",col = "#008080", main = "")

#Correlação de Y com cada Xi
vet_cor_y = cor(dados[,1], dados[,-1])
vet_cor_y
sum(vet_cor_y > 0.9)
sum(vet_cor_y < -0.9)
barplot(sort(vet_cor_y), horiz = TRUE, las = 1,col = ifelse(sort(vet_cor_y) > 0, "#008080", "#FF4500"), border = NA, xlim = c(-1, 1), xlab = "Correlação de cada covariável com Y", main = "")
abline(v = c(max(vet_cor_y), min(vet_cor_y)), lty = 2, col = "gray40")
abline(v = 0, col = "black")
text(x = max(vet_cor_y), y = length(vet_cor_y) + 1, labels = round(max(vet_cor_y), 2), pos = 4, col = "gray40", cex = 0.8)
text(x = min(vet_cor_y), y = length(vet_cor_y) + 1, labels = round(min(vet_cor_y), 2), pos = 2, col = "gray40", cex = 0.8)

#Matriz de Correlação de X
matriz_cor_x = cor(dados[,-1])
det_valor = det(matriz_cor_x)
dci = 1 - det(matriz_cor_x)^(1/(ncol(matriz_cor_x)-1)) #>0.8 a multicolinearidade geral é alta
cor_positivas = sum(matriz_cor_x > 0.9) / 2
cor_negativas = sum(matriz_cor_x < -0.9) / 2
corrplot(matriz_cor_x, method = "color", col = colorRampPalette(c("#FF4500", "white", "#008080"))(200), type = "upper", tl.pos = "n", addCoef.col = NULL, diag = FALSE)
text(x = 0.5, y = 1, labels = paste(sprintf("det = %.3f", det_valor), sprintf("dci = %.3f", dci), sprintf("Corr > 0.9: %.0f", cor_positivas), sprintf("Corr < -0.9: %.0f", cor_negativas), sep = "\n"), pos = 3, cex = 1.2, col = "black", xpd = TRUE)
qgraph(matriz_cor, layout = "spring", edge.color = "#FF4500", node.width = 1.5, label.cex = 1.5, labels = colnames(matriz_cor), label.prop = 0.9, label.scale.equal = TRUE, mar = c(1, 1, 1, 1))
#Separando em conjunto de treino e de teste
treino_ind = sample(seq_len(nrow(dados)), size = floor(0.70*nrow(dados)))
dados_treino = as.data.frame(dados[treino_ind,])
x_treino = as.matrix(cbind(1, dados_treino[,-1]))
y_treino = as.vector(dados_treino[,1])

dados_teste = as.data.frame(dados[-treino_ind,])
x_teste = as.matrix(cbind(1, dados_teste[,-1]))
y_teste = as.vector(dados_teste[,1])

#Modelos
#tempo_lsquares = system.time(modelo_lsquares <- minimos_quadrados(y_treino, x_treino))[3]
#tempo_sma = system.time(modelo_sma <- sma(y_treino, x_treino))[3]
#tempo_backward = system.time(modelo_backward <- backward(y_treino, x_treino))[3]
#tempo_forward = system.time(modelo_forward <- forward(y_treino, x_treino))[3]
#tempo_stepwise = system.time(modelo_stepwise <- stepwise(y_treino, x_treino))[3]
#tempo_lasso = sma_tempo = system.time(modelo_lasso <- lasso(y_treino, x_treino))[3]
#tempo_enet = system.time(modelo_enet <- enet(y_treino, x_treino))[3]
#tempo_gsp_sigmaquaBeta10 = system.time(gsp_sigmaquaBeta10 <- gibbs_sampler(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 10))[3]
#tempo_gsp_sigmaquaBeta50 = system.time(gsp_sigmaquaBeta50 <- gibbs_sampler(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 50))[3]
#tempo_gsp_sigmaquaBeta100 = system.time(gsp_sigmaquaBeta100 <- gibbs_sampler(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 100))[3]
#tempo_gsp_sigmaquaBeta200 = system.time(gsp_sigmaquaBeta200 <- gibbs_sampler(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 200))[3]
#tempo_gsp_sigmaquaBeta500 = system.time(gsp_sigmaquaBeta500 <- gibbs_sampler(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 500))[3]
#tempo_ddrj_sigmaquaBeta10= system.time(ddrj_sigmaquaBeta10 <- ddrj(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 10))[3]
#tempo_ddrj_sigmaquaBeta50= system.time(ddrj_sigmaquaBeta50 <- ddrj(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 50))[3]
#tempo_ddrj_sigmaquaBeta100= system.time(ddrj_sigmaquaBeta100 <- ddrj(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 100))[3]
#tempo_ddrj_sigmaquaBeta200= system.time(ddrj_sigmaquaBeta200 <- ddrj(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 200))[3]
#tempo_ddrj_sigmaquaBeta500= system.time(ddrj_sigmaquaBeta500 <- ddrj(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 500))[3]
#tempo_sslab_sigmaquaBeta10 = system.time(sslab_sigmaquaBeta10 <- spike_slab(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 10))[3]
#tempo_sslab_sigmaquaBeta50= system.time(sslab_sigmaquaBeta50 <- spike_slab(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 50))[3]
#tempo_sslab_sigmaquaBeta100= system.time(sslab_sigmaquaBeta100 <- spike_slab(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 100))[3]
#tempo_sslab_sigmaquaBeta200= system.time(sslab_sigmaquaBeta200 <- spike_slab(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 200))[3]
#tempo_sslab_sigmaquaBeta500= system.time(sslab_sigmaquaBeta500 <- spike_slab(as.vector(dados_treino[,1]), cbind(1, as.matrix(dados_treino[,-1])), 10000, 500))[3]

save.image(file = "/home/giacomini/Testes/ambiente.RData")
load("/home/giacomini/Testes/ambiente.RData")

#Vetores de Betas estimados pelos métodos bayesianos
betas_gsp_sigmaquaBeta10 = rowMeans(gsp_sigmaquaBeta10[[1]])
betas_gsp_sigmaquaBeta10[gsp_sigmaquaBeta10[[3]][,3]==0] = 0
betas_gsp_sigmaquaBeta50 = rowMeans(gsp_sigmaquaBeta50[[1]])
betas_gsp_sigmaquaBeta50[gsp_sigmaquaBeta50[[3]][,3]==0] = 0
betas_gsp_sigmaquaBeta100 = rowMeans(gsp_sigmaquaBeta100[[1]])
betas_gsp_sigmaquaBeta100[gsp_sigmaquaBeta100[[3]][,3]==0] = 0
betas_gsp_sigmaquaBeta200 = rowMeans(gsp_sigmaquaBeta200[[1]])
betas_gsp_sigmaquaBeta200[gsp_sigmaquaBeta200[[3]][,3]==0] = 0
betas_gsp_sigmaquaBeta500 = rowMeans(gsp_sigmaquaBeta500[[1]])
betas_gsp_sigmaquaBeta500[gsp_sigmaquaBeta500[[3]][,3]==0] = 0
betas_sslab_sigmaquaBeta10 = rowMeans(sslab_sigmaquaBeta10[[1]])
betas_sslab_sigmaquaBeta10[rowMeans(sslab_sigmaquaBeta10[[3]])<0.5] = 0
betas_sslab_sigmaquaBeta50 = rowMeans(sslab_sigmaquaBeta50[[1]])
betas_sslab_sigmaquaBeta50[rowMeans(sslab_sigmaquaBeta50[[3]])<0.5] = 0
betas_sslab_sigmaquaBeta100 = rowMeans(sslab_sigmaquaBeta100[[1]])
betas_sslab_sigmaquaBeta100[rowMeans(sslab_sigmaquaBeta100[[3]])<0.5] = 0
betas_sslab_sigmaquaBeta200 = rowMeans(sslab_sigmaquaBeta200[[1]])
betas_sslab_sigmaquaBeta200[rowMeans(sslab_sigmaquaBeta200[[3]])<0.5] = 0
betas_sslab_sigmaquaBeta500 = rowMeans(sslab_sigmaquaBeta500[[1]])
betas_sslab_sigmaquaBeta500[rowMeans(sslab_sigmaquaBeta500[[3]])<0.5] = 0
betas_ddrj_sigmaquaBeta10 = rowMeans(ddrj_sigmaquaBeta10[[1]])
betas_ddrj_sigmaquaBeta10[rowMeans(ddrj_sigmaquaBeta10[[3]])<0.5] = 0
betas_ddrj_sigmaquaBeta50 = rowMeans(ddrj_sigmaquaBeta50[[1]])
betas_ddrj_sigmaquaBeta50[rowMeans(ddrj_sigmaquaBeta50[[3]])<0.5] = 0
betas_ddrj_sigmaquaBeta100 = rowMeans(ddrj_sigmaquaBeta100[[1]])
betas_ddrj_sigmaquaBeta100[rowMeans(ddrj_sigmaquaBeta100[[3]])<0.5] = 0
betas_ddrj_sigmaquaBeta200 = rowMeans(ddrj_sigmaquaBeta200[[1]])
betas_ddrj_sigmaquaBeta200[rowMeans(ddrj_sigmaquaBeta200[[3]])<0.5] = 0
betas_ddrj_sigmaquaBeta500 = rowMeans(ddrj_sigmaquaBeta500[[1]])
betas_ddrj_sigmaquaBeta500[rowMeans(ddrj_sigmaquaBeta500[[3]])<0.5] = 0

indicadora_freq = data.frame(
  "LSquares" = modelo_lsquares[[3]],
  "SMA" = modelo_sma[[3]],
  "StepWise" = modelo_stepwise[[3]],
  "Forward" = modelo_forward[[3]],
  "Backward" = modelo_backward[[3]],
  "LASSO" = modelo_lasso[[3]],
  "ENET" = modelo_enet[[3]]
)

betas_freq = data.frame(
  "LSquares" = modelo_lsquares[[1]],
  "SMA" = modelo_sma[[1]],
  "StepWise" = modelo_stepwise[[1]],
  "Forward" = modelo_forward[[1]],
  "Backward" = modelo_backward[[1]],
  "LASSO" = modelo_lasso[[1]],
  "ENET" = modelo_enet[[1]]
)

indicadora_bayes = data.frame(
  "GSP-10" = gsp_sigmaquaBeta10[[3]][,3],
  "GSP-50" = gsp_sigmaquaBeta50[[3]][,3],
  "GSP-100" = gsp_sigmaquaBeta100[[3]][,3],
  "GSP-200" = gsp_sigmaquaBeta200[[3]][,3],
  "GSP-500" = gsp_sigmaquaBeta500[[3]][,3],
  "SSLAB-10" = (rowMeans(ddrj_sigmaquaBeta10[[3]])>=0.5)*1,
  "SSLAB-50" = (rowMeans(ddrj_sigmaquaBeta50[[3]])>=0.5)*1,
  "SSLAB-100" = (rowMeans(ddrj_sigmaquaBeta100[[3]])>=0.5)*1,
  "SSLAB-200" = (rowMeans(ddrj_sigmaquaBeta200[[3]])>=0.5)*1,
  "SSLAB-500" = (rowMeans(ddrj_sigmaquaBeta500[[3]])>=0.5)*1,
  "DDRJ-10" = (rowMeans(ddrj_sigmaquaBeta10[[3]])>=0.5)*1,
  "DDRJ-50" = (rowMeans(ddrj_sigmaquaBeta50[[3]])>=0.5)*1,
  "DDRJ-100" = (rowMeans(ddrj_sigmaquaBeta100[[3]])>=0.5)*1,
  "DDRJ-200" = (rowMeans(ddrj_sigmaquaBeta200[[3]])>=0.5)*1,
  "DDRJ-500" = (rowMeans(ddrj_sigmaquaBeta500[[3]])>=0.5)*1
)

betas_bayes = data.frame(
  "GSP-10" = betas_gsp_sigmaquaBeta10,
  "GSP-50" = betas_gsp_sigmaquaBeta50,
  "GSP-100" = betas_gsp_sigmaquaBeta100,
  "GSP-200" = betas_gsp_sigmaquaBeta200,
  "GSP-500" = betas_gsp_sigmaquaBeta500,
  "SSLAB-10" = betas_sslab_sigmaquaBeta10,
  "SSLAB-50" = betas_sslab_sigmaquaBeta50,
  "SSLAB-100" = betas_sslab_sigmaquaBeta100,
  "SSLAB-200" = betas_sslab_sigmaquaBeta200,
  "SSLAB-500" = betas_sslab_sigmaquaBeta500,
  "DDRJ-10" = betas_ddrj_sigmaquaBeta10,
  "DDRJ-50" = betas_ddrj_sigmaquaBeta50,
  "DDRJ-100" = betas_ddrj_sigmaquaBeta100,
  "DDRJ-200" = betas_ddrj_sigmaquaBeta200,
  "DDRJ-500" = betas_ddrj_sigmaquaBeta500
)

metricas = data.frame(
  "LSquares" = c(modelo_lsquares[[2]], sum((y_teste-x_teste%*%modelo_lsquares[[1]])^2)/length(y_teste), sum(abs(y_teste-x_teste%*%modelo_lsquares[[1]]))/length(y_teste), tempo_lsquares),
  "SMA" = c(modelo_sma[[2]], sum((y_teste-x_teste%*%modelo_sma[[1]])^2)/length(y_teste), sum(abs(y_teste-x_teste%*%modelo_sma[[1]]))/length(y_teste), tempo_sma),
  "StepWise" = c(modelo_stepwise[[2]], sum((y_teste-x_teste%*%modelo_stepwise[[1]])^2)/length(y_teste), sum(abs(y_teste-x_teste%*%modelo_stepwise[[1]]))/length(y_teste), tempo_stepwise),
  "Forward" = c(modelo_forward[[2]], sum((y_teste-x_teste%*%modelo_forward[[1]])^2)/length(y_teste), sum(abs(y_teste-x_teste%*%modelo_forward[[1]]))/length(y_teste), tempo_forward),
  "Backward" = c(modelo_backward[[2]], sum((y_teste-x_teste%*%modelo_backward[[1]])^2)/length(y_teste), sum(abs(y_teste-x_teste%*%modelo_backward[[1]]))/length(y_teste), tempo_backward),
  "LASSO" = c(modelo_lasso[[2]], sum((y_teste-x_teste%*%modelo_lasso[[1]])^2)/length(y_teste), sum(abs(y_teste-x_teste%*%modelo_lasso[[1]]))/length(y_teste), tempo_lasso),
  "ENET" = c(modelo_enet[[2]], sum((y_teste-x_teste%*%modelo_enet[[1]])^2)/length(y_teste), sum(abs(y_teste-x_teste%*%modelo_enet[[1]]))/length(y_teste), tempo_enet),
  "GSP-10" = c(mean(gsp_sigmaquaBeta10[[2]]), sum((y_teste-x_teste%*%betas_gsp_sigmaquaBeta10)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_gsp_sigmaquaBeta10))/length(y_teste), tempo_gsp_sigmaquaBeta10),
  "GSP-50" = c(mean(gsp_sigmaquaBeta50[[2]]), sum((y_teste-x_teste%*%betas_gsp_sigmaquaBeta50)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_gsp_sigmaquaBeta50))/length(y_teste), tempo_gsp_sigmaquaBeta50),
  "GSP-100" = c(mean(gsp_sigmaquaBeta100[[2]]), sum((y_teste-x_teste%*%betas_gsp_sigmaquaBeta100)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_gsp_sigmaquaBeta100))/length(y_teste), tempo_gsp_sigmaquaBeta100),
  "GSP-200" = c(mean(gsp_sigmaquaBeta200[[2]]), sum((y_teste-x_teste%*%betas_gsp_sigmaquaBeta200)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_gsp_sigmaquaBeta200))/length(y_teste), tempo_gsp_sigmaquaBeta200),
  "GSP-500" = c(mean(gsp_sigmaquaBeta500[[2]]), sum((y_teste-x_teste%*%betas_gsp_sigmaquaBeta500)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_gsp_sigmaquaBeta500))/length(y_teste), tempo_gsp_sigmaquaBeta500),
  "SSLAB-10" = c(mean(sslab_sigmaquaBeta10[[2]]), sum((y_teste-x_teste%*%betas_sslab_sigmaquaBeta10)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_sslab_sigmaquaBeta10))/length(y_teste), tempo_sslab_sigmaquaBeta10),
  "SSLAB-50" = c(mean(sslab_sigmaquaBeta50[[2]]), sum((y_teste-x_teste%*%betas_sslab_sigmaquaBeta50)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_sslab_sigmaquaBeta50))/length(y_teste), tempo_sslab_sigmaquaBeta50),
  "SSLAB-100" = c(mean(sslab_sigmaquaBeta100[[2]]), sum((y_teste-x_teste%*%betas_sslab_sigmaquaBeta100)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_sslab_sigmaquaBeta100))/length(y_teste), tempo_sslab_sigmaquaBeta100),
  "SSLAB-200" = c(mean(sslab_sigmaquaBeta200[[2]]), sum((y_teste-x_teste%*%betas_sslab_sigmaquaBeta200)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_sslab_sigmaquaBeta200))/length(y_teste), tempo_sslab_sigmaquaBeta200),
  "SSLAB-500" = c(mean(sslab_sigmaquaBeta500[[2]]), sum((y_teste-x_teste%*%betas_sslab_sigmaquaBeta500)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_sslab_sigmaquaBeta500))/length(y_teste), tempo_sslab_sigmaquaBeta500),
  "DDRJ-10" = c(mean(ddrj_sigmaquaBeta10[[2]]), sum((y_teste-x_teste%*%betas_ddrj_sigmaquaBeta10)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_ddrj_sigmaquaBeta10))/length(y_teste), tempo_ddrj_sigmaquaBeta10),
  "DDRJ-50" = c(mean(ddrj_sigmaquaBeta50[[2]]), sum((y_teste-x_teste%*%betas_ddrj_sigmaquaBeta50)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_ddrj_sigmaquaBeta50))/length(y_teste), tempo_ddrj_sigmaquaBeta50),
  "DDRJ-100" = c(mean(ddrj_sigmaquaBeta100[[2]]), sum((y_teste-x_teste%*%betas_ddrj_sigmaquaBeta100)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_ddrj_sigmaquaBeta100))/length(y_teste), tempo_ddrj_sigmaquaBeta100),
  "DDRJ-200" = c(mean(ddrj_sigmaquaBeta200[[2]]), sum((y_teste-x_teste%*%betas_ddrj_sigmaquaBeta200)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_ddrj_sigmaquaBeta200))/length(y_teste), tempo_ddrj_sigmaquaBeta200),
  "DDRJ-500" = c(mean(ddrj_sigmaquaBeta500[[2]]), sum((y_teste-x_teste%*%betas_ddrj_sigmaquaBeta500)^2)/length(y_teste), sum(abs(y_teste-x_teste%*%betas_ddrj_sigmaquaBeta500))/length(y_teste), tempo_ddrj_sigmaquaBeta500)
)
rownames(metricas) = c("EQMTREINO", "EQMTESTE", "EAMTESTE", "TEMPO")
metricas = cbind(t(metricas), c(colSums(indicadora_freq), colSums(indicadora_bayes)))
colnames(metricas)[5] = "P-SELECT"

arquivo = createWorkbook()
nome_arquivo = "dados_reais.xlsx"

addWorksheet(arquivo, "indicadora_freq")
addWorksheet(arquivo, "betas_freq")
addWorksheet(arquivo, "indicadora_bayes")
addWorksheet(arquivo, "betas_bayes")
addWorksheet(arquivo, "metricas")
writeData(arquivo, sheet = "indicadora_freq", indicadora_freq, rowNames = F, colNames = T)
writeData(arquivo, sheet = "betas_freq", betas_freq, rowNames = F, colNames = T)
writeData(arquivo, sheet = "indicadora_bayes", indicadora_bayes, rowNames = F, colNames = T)
writeData(arquivo, sheet = "betas_bayes", betas_bayes, rowNames = F, colNames = T)
writeData(arquivo, sheet = "metricas", metricas, rowNames = T, colNames = T)
saveWorkbook(arquivo, file.path("/home/giacomini/Testes", nome_arquivo), overwrite = TRUE)