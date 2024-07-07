adjust_model_env <- function(mod, env, script.dir) {
  if("highH2" %in% env) {
    env_dt <- fread(paste0(script.dir, "/../dat/env/env_highH2.tsv"), header = F)
    env_dt[, V1 := paste0(V1,"_c0")]
    env_dt <- env_dt[V1 %in% mod@react_id]
    if(nrow(env_dt) > 0) {
      for(i in 1:nrow(env_dt)) {
        if(env_dt[i, V2] == ">")
          mod <- changeBounds(mod, react = env_dt[i, V1], lb = 0, ub = COBRAR_SETTINGS("MAXIMUM"))
        if(env_dt[i, V2] == "<")
          mod <- changeBounds(mod, react = env_dt[i, V1], lb = -COBRAR_SETTINGS("MAXIMUM"), ub = 0)
        if(env_dt[i, V2] == "=")
          mod <- changeBounds(mod, react = env_dt[i, V1], lb = -COBRAR_SETTINGS("MAXIMUM"), ub = COBRAR_SETTINGS("MAXIMUM"))
      }
    }
  }
  
  return(mod)
}