# jobs to run over weekend
setwd(bootResWd)
cur_num  = 350
cur_cpus = 7

fun_extr = function(cur_mod) {
  ## extraction function to return fixeffect and mse
  cur_fe  = fixef(cur_mod)
  cur_mse = mean((1/(1+exp(-predict(cur_mod)))-as.numeric(as.character(cur_mod@frame$accept_reject)))^2)
  cur_fe  = c(cur_fe,cur_mse)
  names(cur_fe)[length(cur_fe)] = 'mse'
  return(cur_fe)
}

# PERMUTATION TESTS
# bootstrap p-value modla_0g (permutation)
effects_under_0_0g = agk.boot.p.mermod(mermod = modla_0g,mermod0 = modla_01,num_cpus = cur_cpus,num = cur_num,fun_extract = fun_extr,cur_control = cur_control,permvars = c('HCPG'),type='perm')
save(file= 'effects_under_0_0g_perm_1000_wc.RData',list=c('effects_under_0_0g'))

# bootstrap p-value modla_cg (permutation)
effects_under_0_cg = agk.boot.p.mermod(mermod = modla_cg,mermod0 = modla_01,num_cpus = cur_cpus,num = cur_num,fun_extract = fun_extr,cur_control = cur_control,permvars = c('cat'),type='perm')
save(file = 'effects_under_0_cg_perm_1000_wc.RData',list=c('effects_under_0_cg'))

# CFINTS
# bootstrap cfint modla_0g (np boot)
boot_cfint_0g = agk.boot.cfint.mermod(mermod = modla_0g,num_cpus = cur_cpus,num = cur_num,fun_extract = fixef,cur_control = cur_control,type = 'non-parametric')
save(file = 'boot_cfint_0g_1000_wc.RData',list=c('boot_cfint_0g'))

# bootstrap cfint modla_cg (np boot)
boot_cfint_cg = agk.boot.cfint.mermod(mermod = modla_cg,num_cpus = cur_cpus,num = cur_num,fun_extract = fixef,cur_control = cur_control,type = 'non-parametric')
save(file = 'boot_cfint_cg_1000_wc.RData',list=c('boot_cfint_cg'))




# CFINT (cgi)
# bootstrap cfint modla_cgi (np boot)
boot_cfint_cgi = agk.boot.cfint.mermod(mermod = modla_cgi,num_cpus = cur_cpus,num = cur_num,fun_extract = fixef,cur_control = cur_control,type = 'non-parametric')
save(file = 'boot_cfint_cgi_1000_wc.RData',list=c('boot_cfint_cgi'))

# PERMUTATION TEST (cgi)
# bootstrap p-value modla_cgi (permutation) (actually not to do cause no sig model comparison anyways)
cur_permvars        = names(fixef(modla_cgirf))
cur_permvars        = cur_permvars[grep('.cat',cur_permvars,fixed=T)]
effects_under_0_cgi = agk.boot.p.mermod(mermod = modla_cgirf,mermod0 = modla_01,num_cpus = cur_cpus,num = cur_num,fun_extract = fun_extr,cur_control = cur_control,permvars = cur_permvars,type='perm')
save(file = 'effects_under_0_cgi_perm_1000_wc.RData',list=c('effects_under_0_cgi'))
