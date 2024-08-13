get_scenario <- function(Risk, Amount, Params, t) {
  
  # get actual proportions cleared (annualised)
  # divde by 11 to annualise since the risk estimate is for 11 years of clearing
  p <- (ifelse(Risk == 0, 0.000000000000001, ifelse(Risk == 1, 0.999999999999999, Risk)) / 11) %>% as.matrix()
  
  # get the baseline
  b <- Amount %>% as.matrix()
  
  get_area_effect <- function(p, b, t, Params) {
    # get the predicted area effects
    Out <- (((((exp(log(p / (1 - p))) / (1 + exp(log(p / (1 - p)))))) * b) - (((exp(log(p / (1 - p)) - Params[1] - Params[2] * t) / (1 + exp(log(p / (1 - p)) - Params[1] - Params[2] * t)))) * b)) * 25 * 25 / 10000) %>% as_tibble()
    OutSums <- Out %>% summarise(across(everything(), \(x) sum(x, na.rm = TRUE)))
    return(OutSums)
  }
  
  # get the output
  Output <- apply(Params, MARGIN = 1, FUN = function(x) {get_area_effect(p, b, t, x)})
  Output <- do.call("rbind", Output)
  return(Output)
}
