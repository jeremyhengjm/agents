#' @title Transition probability of SIS agent-based model
#' @description Computes transition probability of SIS agent-based model.
#' @param model_config a list containing: 
#' \itemize{
#' \item 'lambda': vector of infection rates of length N 
#' \item 'gamma': vector of recovery rates of length N 
#' \item 'network_type': character "full" or "net" specifying if network is fully connected or has a given network structure
#' \item 'adjacency': matrix indicating interaction of agents of dimension N x N (if model_config$network_type == "net)
#' \item 'degree' vector of node degree of length N (if model_config$network_type == "net)
#' }
#' @param X matrix of agent states of dimension N x P 
#' @param I vector infected agents of length P (to avoid additional computation if this summary is needed)
#' @return matrix of transition probabilities of dimension N x P 
#' @export  

sis_alpha_prob <- function(model_config, X, I){
  # fully connected network
  if (model_config$network_type == "full"){
    alpha <- sis_alpha_full(X, I, model_config$lambda, model_config$gamma)
  } 
  
  # network specified by adjacency matrix and degree vector
  if (model_config$network_type == "net"){
    alpha <- sis_alpha_net(X, model_config$lambda, model_config$gamma, model_config$adjacency, model_config$degree)
  }
  
  return(alpha)
}