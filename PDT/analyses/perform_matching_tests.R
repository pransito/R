# two functions to perform the matching test
# first the dat_match has to be split and packed into 
# Cohort data frames; secondly the matching is performed on those
agk.perform.matching.splitdfs = function(dat_match,cur_study = NULL) {
  # what are the names of the to-be-matched vars?
  cur_names = c('edu_years','edu_years_voca','edu_hollingshead', 'income_personal','debt_overall_num','debt_gambling_num','smoking_ftdt','Age',               
                'audit','BDI','SOGS','KFG','BIS','GBQ_persi','GBQ_illus','BIG')        
  
  # splitting
  dat_match_MRI              = subset(dat_match,Cohort == "MRT")
  dat_match_POSTPILOT        = subset(dat_match,(Cohort == "PGPilot" | Cohort == "POSTPILOT"))
  dat_match_PP_MRI           = subset(dat_match,(Cohort == "PGPilot" | Cohort == "POSTPILOT" | Cohort == "MRT"))
  dat_match_PP_MRI$study     = NA
  dat_match_PP_MRI$study[dat_match_PP_MRI$Cohort == "PGPilot"] = "POSTPILOT"
  dat_match_PP_MRI$study[dat_match_PP_MRI$Cohort == "POSTPILOT"] = "POSTPILOT"
  dat_match_PP_MRI$study[dat_match_PP_MRI$Cohort == "MRT"] = "MRT"
  dat_match_PP_MRI$study = as.factor(dat_match_PP_MRI$study)
  
  # packing
  dfs          = list(dat_match_MRI,dat_match_POSTPILOT,dat_match_POSTPILOT,dat_match_PP_MRI)
  dfs[[3]]     = subset(dfs[[3]],HCPG == "PG") # for gender study only use PG subjects
  cur_groups   = c("HCPG","HCPG","dem_gender","study")
  cur_matching = c("MRT","POSTPILOT","POSTPILOT_GENDER","PP_MRI")
  
  if (!is.null(cur_study)) {
    # only if current study is given
    # select only the Cohort of interest
    if (cur_study == 'POSTPILOT') {
      cur_ind =2
    } else if (cur_study == 'MRT') {
      cur_ind = 1
    } else {
      stop(paste('rec_matching_test not implemented for',cur_study,'study.'))
    }
    dfs          = dfs[cur_ind]
    cur_groups   = cur_groups[cur_ind]
    cur_matching = cur_matching[cur_ind]
  }
  
  return(list(dfs=dfs,cur_groups=cur_groups,cur_matching=cur_matching,cur_names = cur_names))
}

agk.perform.matching.tests = function(dfs,cur_groups,cur_matching,path_mtc,write_match_tables,cur_names) {
  # running the matching tests based on dat_match's
  # packed in df's
  # create a list where the resulting matching tables are stored
  match_result_tables = list()
  
  for (kk in 1:length(cur_groups)) {
    df        = dfs[[kk]]
    cur_group = cur_groups[kk]
    cur_test  = list()
    
    for (ii in 1:length(cur_names)) {
      cur_test[[ii]] = f.difftest(df,cur_names[ii],cur_group)
    }
    
    cur_match_test = cur_test[[1]]
    for (ii in 2:length(cur_test)) {
      cur_match_test = rbind(cur_match_test,cur_test[[ii]])
    }
    cur_match_test$cur_var = as.character(cur_match_test$cur_var)
    cur_match_test$matched = as.character(cur_match_test$matched)
    if (match_boot == 1) {
      cur_match_test$t_matched_boot = as.character(cur_match_test$t_matched_boot)
    } else if (match_boot == 2) {
      cur_match_test$t_matched_perm = as.character(cur_match_test$t_matched_perm)
    } else {
      cur_match_test$t_matched = as.character(cur_match_test$t_matched)
    }
    
    
    # Gender
    if (cur_group != "dem_gender") {
      cur_formula = as.formula(paste0("~",cur_group,"+dem_gender"))
      tmp = xtabs(cur_formula,data = df)
      tmp = tmp[c(1,2),c(1,2)]
      rfh = as.numeric(round(tmp[rownames(tmp) == rownames(tmp)[1],colnames(tmp) == "weiblich"]/(tmp[1,1]+tmp[1,2]),2))
      rfg = as.numeric(round(tmp[rownames(tmp) == rownames(tmp)[2],colnames(tmp) == "weiblich"]/(tmp[2,1]+tmp[2,2]),2))
      tmp = chisq.test(tmp,simulate.p.value = T)
      cur_matched = ifelse(tmp$p.value < m_crit,"NO","YES")
      gline = c("ratio_female",rfh,NA,rfg,NA,NA,tmp$p.value,cur_matched,NA,NA,cur_matched)
      cur_match_test = rbind(cur_match_test,gline)
    }
    
    # Occupation
    cur_formula = as.formula(paste0("~",cur_group,"+unemployed"))
    tmp = xtabs(cur_formula,data = df)
    rfh = as.numeric(round(tmp[rownames(tmp) == rownames(tmp)[1],colnames(tmp) == "unemployed"]/(tmp[1,1]+tmp[1,2]),2))
    rfg = as.numeric(round(tmp[rownames(tmp) == rownames(tmp)[2],colnames(tmp) == "unemployed"]/(tmp[2,1]+tmp[2,2]),2))
    tmp = chisq.test(tmp,simulate.p.value = T)
    cur_matched = ifelse(tmp$p.value < m_crit,"NO","YES")
    gline = c("ratio_unemployed",rfh,NA,rfg,NA,NA,tmp$p.value,cur_matched,NA,NA,cur_matched)
    cur_match_test = rbind(cur_match_test,gline)
    
    # Smoking
    cur_formula = as.formula(paste0("~",cur_group,"+smoking"))
    tmp = xtabs(cur_formula,data = df)
    rfh = round(tmp[rownames(tmp) == rownames(tmp)[1],colnames(tmp) == "ja"]/(tmp[1,1]+tmp[1,2]),2)
    rfg = round(tmp[rownames(tmp) == rownames(tmp)[2],colnames(tmp) == "ja"]/(tmp[2,1]+tmp[2,2]),2)
    tmp = chisq.test(tmp,simulate.p.value = T)
    cur_matched = ifelse(tmp$p.value < m_crit,"NO","YES")
    gline = c("smoking_yes",rfh,NA,rfg,NA,NA,tmp$p.value,cur_matched,NA,NA,cur_matched)
    cur_match_test = rbind(cur_match_test,gline)
    
    # handedness
    cur_formula = as.formula(paste0("~",cur_group,"+handedness"))
    tmp  = xtabs(cur_formula,data = df,na.action = NULL)
    tmp  = tmp[c(1:2),c(1:2)]
    rfh = round(tmp[rownames(tmp) == rownames(tmp)[1],colnames(tmp) == "rechts"]/(sum(tmp[1,])),2)
    rfg = round(tmp[rownames(tmp) == rownames(tmp)[2],colnames(tmp) == "rechts"]/(sum(tmp[2,])),2)
    tmp = chisq.test(tmp,simulate.p.value = T)
    cur_matched = ifelse(tmp$p.value < m_crit,"NO","YES")
    gline = c("ratio_righthandedness",rfh,NA,rfg,NA,NA,tmp$p.value,cur_matched,NA,NA,cur_matched)
    cur_match_test = rbind(cur_match_test,gline)
    
    # add N
    cur_formula = as.formula(paste0("~",cur_group))
    tmp = xtabs(cur_formula,data = df)
    rfh = tmp[1]
    rfg = tmp[2]
    cur_matched = ifelse(abs((rfh - rfg)) > 1,"NO","YES")
    gline = c("N",rfh,NA,rfg,NA,NA,NA,cur_matched,NA,NA,cur_matched)
    cur_match_test = rbind(cur_match_test,gline)
    
    # add info, of how many good physio data is there
    # getting the intersection of people with data and with physio data
    df_phys = df[df[[which(names(df) =="VPPG")]] %in% subs_good_phys,]
    
    # add N_phys
    cur_formula = as.formula(paste0("~",cur_group))
    tmp = xtabs(cur_formula,data = df_phys)
    rfh = tmp[1]
    rfg = tmp[2]
    cur_matched = ifelse(abs((rfh - rfg)) > 1,"NO","YES")
    gline = c("N_phys",rfh,NA,rfg,NA,NA,NA,cur_matched,NA,NA,cur_matched)
    cur_match_test = rbind(cur_match_test,gline)
    
    # add the needs matching info
    if (cur_group != "dem_gender") {
      needs_matching                      = c("YES","YES","YES","YES","DM","NO","YES","YES","YES","DM","NO","NO","DM","DM","DM",'NO',"YES","YES","YES","YES","DM","DM")
    } else {
      needs_matching                      = c("YES","YES","YES","YES","YES","NO","YES","YES","YES","YES","YES","YES","YES","DM","DM",'NO',"YES","YES","YES","DM","DM")
    }
    
    cur_match_test$needs_matching       = needs_matching
    cur_match_test$ok_or_checkneeded    = NA
    for (ll in 1:length(cur_match_test$needs_matching)) {
      if (needs_matching[ll] == "DM") {
        cur_match_test$ok_or_checkneeded[ll] = "DOESNT_MATTER"
      } else if (cur_match_test$needs_matching[ll] == cur_match_test$matched[ll]) {
        cur_match_test$ok_or_checkneeded[ll] = "OKAY"
      } else {
        cur_match_test$ok_or_checkneeded[ll] = "NEEDS_CHECK"
      }
    }
    
    # rounding for nicer display
    cur_p_vars = names(cur_match_test)[grep('p-value',names(cur_match_test))]
    cur_r_vars = c(cur_p_vars,'cur_se')
    for (gg in 1:length(cur_p_vars)) {
      cur_match_test[[cur_p_vars[gg]]] = round(as.numeric(cur_match_test[[cur_p_vars[gg]]]),digits = 3)
    }
    
    # store
    match_result_tables[[kk]] = cur_match_test
    
    if (write_match_tables) {
      write.xlsx(cur_match_test, paste0(path_mtc,"/",cur_matching[kk],"_matching.xlsx"),row.names=FALSE, showNA=TRUE)
    }
    
  }
  # return
  return(match_result_tables)
}

agk.interpolating.dat_match = function(dfs,cur_groups,cur_names,cur_gr_levs) {
  # function to interpolate using mean
  for (ii in 1:length(dfs)) {
    # loop over studies (dfs)
    cur_df = dfs[[ii]]
    # get group specific dfs
    cur_df_A = cur_df[cur_df[[cur_groups[[ii]]]] == cur_gr_levs[[ii]][1],]
    cur_df_B = cur_df[cur_df[[cur_groups[[ii]]]] == cur_gr_levs[[ii]][2],]
    
    for (jj in 1:length(cur_names)) {
      cur_df_A[[cur_names[jj]]][is.na(cur_df_A[[cur_names[jj]]])] = mean.rmna(cur_df_A[[cur_names[jj]]])
      cur_df_B[[cur_names[jj]]][is.na(cur_df_B[[cur_names[jj]]])] = mean.rmna(cur_df_B[[cur_names[jj]]])
    }
    # pack pack into cur_df/dfs
    cur_df[cur_df[[cur_groups[[ii]]]] == cur_gr_levs[[ii]][1],] = cur_df_A
    cur_df[cur_df[[cur_groups[[ii]]]] == cur_gr_levs[[ii]][2],] = cur_df_B
    dfs[[ii]] = cur_df
  }
  return(dfs)
}

agk.domatch = function(which_studies,desired_n,dfs,cur_groups,cur_names_dom) {
  # function to pick as many as n_desired subjects in PG of which_studies cohorts
  # to find the closest HC subject according to cur_names variables
  # if desired_n is lower than available subjects than it will drop subjects in the end
  # will take all PGs: find the optimal HC partners and then drop the worst couples
  # until desired_n is reached
  # ORDER of PGs IS DISREGARDED HERE (NOT FEASABLE)
  # however note: if fewer PGs than HCs then goes from PGs, if vice versa then goes from
  # HCs; always does all couples starting from the smaller group, then picks best couples

  # normalizes data between 0 and 1 beforehand
  
  # drop lists
  drop_lists = list()
  for(ii in 1:length(which_studies)) {
    # loop over studies
    cur_df = dfs[[ii]]
    
    # transform desired vars to numeric and between 0 and 1
    cur_df_trnsf = cur_df
    cur_df_trnsf[cur_names_dom] = as.data.frame(lapply(cur_df_trnsf[cur_names_dom],FUN=as.numeric))
    cur_df_trnsf[cur_names_dom] = as.data.frame(lapply(cur_df_trnsf[cur_names_dom],FUN=scale))
    
    # weight handedness in the MRT study
    if (which_studies[ii] =='MRT') {
      cur_df_trnsf$handedness       = cur_df_trnsf$handedness * 2
      cur_df_trnsf$Age              = cur_df_trnsf$Age * 4
      cur_df_trnsf$edu_hollingshead = cur_df_trnsf$edu_hollingshead * 1.5
    }

    # get group specific dfs
    cur_df_HC      = cur_df_trnsf[cur_df_trnsf[[cur_groups[[ii]]]] == 'HC',]
    cur_df_PG      = cur_df_trnsf[cur_df_trnsf[[cur_groups[[ii]]]] == 'PG',]
    cur_df_HC_bcp  = cur_df_HC
    cur_df_PG_bcp  = cur_df_PG
    
    # check numbers
    if (desired_n[[ii]][2] > length(cur_df_PG[,1])) {
      # PG desired_n bigger than available n of PG 
      disp("Not as many PGs available as desired.")
      next
    } else if (desired_n[[ii]][2] == length(cur_df_PG[,1])) {
      # PG desired_n equal to available n of PG
      if (desired_n[[ii]][1] >= length(cur_df_HC[,1])) {
        # HC is equal or smaller than desired n
        disp('Person-wise matching not possible because too few HC subjects to choose from.')
        next
      }
    } else if (desired_n[[ii]][2] < length(cur_df_PG[,1])) {
      # desired n of PG is smaller than available PGs
      # This case is always ok (if there are HCs)
    }
    
    # check from which group to start
    if (length(cur_df_HC_bcp[,1]) < length(cur_df_PG_bcp[,1])) {
      switch_groups  = T
      cur_df_HCb     = cur_df_HC
      cur_df_HC_bcpb = cur_df_HC_bcp
      cur_df_HC      = cur_df_PG
      cur_df_HC_bcp  = cur_df_PG_bcp
      cur_df_PG      = cur_df_HCb
      cur_df_PG_bcp  = cur_df_HC_bcpb
    } else if (length(cur_df_HC_bcp[,1]) > length(cur_df_PG_bcp[,1])) {
      switch_groups = F
    } else {
      switch_groups = F
    }
    
    # get the desired vars
    cur_df_HC_vars = cur_df_HC[cur_names_dom]
    
    # transform to numeric [REDUNDANT?]
    cur_df_HC_vars = as.data.frame(lapply(cur_df_HC_vars,FUN=as.numeric))
    
    agk.domatch.get.distances = function(cur_df_PG,cur_df_HC,cur_df_HC_vars,desired_n,ii) {
      # prepare to collect the optimal distances per coupling
      cur_dists  = c()
      chosen_HCs = c()
      chosen_PGs = c()
      for (jj in 1:length(cur_df_PG[,1])) {
        # loop over PG subjects (or HC subjects whoever is fewer;)
        # get this subject's matching parameter vector
        cur_sub = cur_df_PG[jj,]
        cur_sub = cur_sub[cur_names_dom]
        # transform to numeric
        cur_sub = as.data.frame(lapply(cur_sub,FUN=as.numeric))
        # make the matrix with first line being the HC sub, all following the PGs
        cur_mat   = rbind(as.matrix(cur_sub),as.matrix(cur_df_HC_vars))
        cur_dist  = as.matrix(dist(cur_mat,method = "euclidean"))
        cur_dist  = cur_dist[-1,1]
        match_sub = which((min(as.numeric(cur_dist)) == as.numeric(cur_dist)))[1]
        # note down the subject chosen and the distance
        # print(min(as.numeric(cur_dist)))
        cur_dists[jj]  = min(as.numeric(cur_dist))
        chosen_PGs[jj] = cur_df_PG$VPPG[jj]
        chosen_HCs[jj] = cur_df_HC$VPPG[match_sub]
        # delete the subject from possible PGs
        cur_df_HC      = cur_df_HC[-match_sub,]
        # preserve data frame structure needed if only one variable
        cur_df_HC_vars = as.data.frame(cur_df_HC_vars[-match_sub,])
        names(cur_df_HC_vars) = cur_names_dom
      }
      
      # take the best n couples
      cur_n = desired_n[[ii]][2]
      chosen_HCsPGs = data.frame(chosen_HCs,chosen_PGs)
      while (length(chosen_HCsPGs[,1]) > cur_n) {
        cur_max       = max(cur_dists)
        cur_ind       = which(cur_max == cur_dists)[1]
        chosen_HCsPGs = chosen_HCsPGs[-cur_ind,]
        cur_dists     = cur_dists[-cur_ind]
      }
      
      # note down which PGs, HCs, were not chosen
      dropped_HCs_m = cur_df_HC_bcp$VPPG[!cur_df_HC_bcp$VPPG %in% chosen_HCsPGs$chosen_HCs]
      dropped_PGs_m = cur_df_PG_bcp$VPPG[!cur_df_PG_bcp$VPPG %in% chosen_HCsPGs$chosen_PGs]
      return(list(dropped_HCs_m = dropped_HCs_m, dropped_PGs_m = dropped_PGs_m,
                  chosen_HCs = chosen_HCsPGs$chosen_HCs, chosen_PGs = chosen_HCsPGs$chosen_PGs,
                  mean_distance = mean(cur_dists)))
    }
    
    # get the coupling up to desired n
    res = agk.domatch.get.distances(cur_df_PG,cur_df_HC,cur_df_HC_vars,desired_n,ii)
    
    # unpack (and invert)
    if (switch_groups) {
      dropped_PGs_m = res$dropped_HCs_m
      dropped_HCs_m = res$dropped_PGs_m
      chosen_PGs    = res$chosen_HCs
      chosen_HCs    = res$chosen_PGs
    } else {
      dropped_HCs_m = res$dropped_HCs_m
      dropped_PGs_m = res$dropped_PGs_m
      chosen_HCs    = res$chosen_HCs
      chosen_PGs    = res$chosen_PGs
    }
    
    # saving the dropped PGs and HCs
    cur_drop_list = list()
    cur_drop_list$dropped_HC_matching = dropped_HCs_m
    cur_drop_list$dropped_PG_matching = dropped_PGs_m
    drop_lists[[ii]]                  = cur_drop_list
    
    # packing it back
    if (switch_groups) {
      cur_df_HC = cur_df_PG_bcp[cur_df_PG_bcp$VPPG %in% chosen_HCs,]
      cur_df_PG = cur_df_HC_bcp[cur_df_HC_bcp$VPPG %in% chosen_PGs,]
    } else {
      cur_df_PG = cur_df_PG_bcp[cur_df_PG_bcp$VPPG %in% chosen_PGs,]
      cur_df_HC = cur_df_HC_bcp[cur_df_HC_bcp$VPPG %in% chosen_HCs,]
    }
    
    # getting the untransformed cur_df
    cur_df_HCut = cur_df[cur_df[[cur_groups[[ii]]]] == 'HC',]
    cur_df_PGut = cur_df[cur_df[[cur_groups[[ii]]]] == 'PG',]
    
    # choosing the chosen PGs and HCs
    cur_df_HCut = cur_df_HCut[cur_df_HCut$VPPG %in% cur_df_HC$VPPG,]
    cur_df_PGut = cur_df_PGut[cur_df_PGut$VPPG %in% cur_df_PG$VPPG,]
    cur_df      = rbind(cur_df_HCut,cur_df_PGut)
    dfs[[ii]] = cur_df
  }
  return(list(dfs = dfs, dropped_HCs_PGs = drop_lists))
}

agk.domatch.elim = function(which_studies,desired_n,dfs,cur_groups,cur_names_dom_narrowed) {
  # DOES NOT WORK WELL DO NOT USE!!!
  # function that takes a dat_match after person-by person matching
  # will try to improve matching further by cutting out couples (HC,PG)
  # which have large distance
  # then tests if matching ok
  # if not continues. else returns
  # assumes equal group sizes; otherwise does nothing (then do agk.domatch first)
  
  dfs_bcp    = dfs
  drop_lists = list()
  for(ii in 1:length(which_studies)) {
    # loop over studies
    cur_df = dfs[[ii]]
    
    # get the available subjects
    h_PG = sum(cur_df$HCPG == 'PG')
    h_HC = sum(cur_df$HCPG == 'HC')
    
    # only if equal group size
    if (h_PG != h_HC) {
      return(...)
    }
    
    # transform desired vars to numeric and scale
    cur_df_trnsf = cur_df
    cur_df_trnsf[cur_names_dom_narrowed] = as.data.frame(lapply(cur_df_trnsf[cur_names_dom_narrowed],FUN=as.numeric))
    cur_df_trnsf[cur_names_dom_narrowed] = as.data.frame(lapply(cur_df_trnsf[cur_names_dom_narrowed],FUN=scale))
    
    # get group specific dfs
    cur_df_HC      = cur_df_trnsf[cur_df_trnsf[[cur_groups[[ii]]]] == 'HC',]
    cur_df_PG      = cur_df_trnsf[cur_df_trnsf[[cur_groups[[ii]]]] == 'PG',]
    cur_df_HC_bcp  = cur_df_HC
    cur_df_PG_bcp  = cur_df_PG
    switch_groups = F
    
    # get the desired vars
    cur_df_HC_vars = cur_df_HC[cur_names_dom_narrowed]
    cur_df_PG_vars = cur_df_PG[cur_names_dom_narrowed]
    
    # get the distances
    # CAREFUL late singles are always in disadvantage because market of prospective partner is empty
    # that is why we order first: singles which have low distance in general come last; difficult ones first
    # using mean distance (works indeed better; mean distance in the end smaller than when doing easiest first)
    HC_kill    = c()
    PG_kill    = c()
    ct         = 0
    # check first if this is even needed
    res        = agk.perform.matching.tests(dfs,cur_groups,cur_matching,path_mtc,0,cur_names)
    score      = res[[ii]]$ok_or_checkneeded[res[[ii]]$cur_var %in% cur_names_dom_narrowed]
    contn_elim = any(score == 'NEEDS_CHECK')
    
    while (contn_elim) {
      disp(paste('Cutting couple number',ct+1))
      cur_dist = distmat(as.matrix(cur_df_HC_vars), as.matrix(cur_df_PG_vars))
      cur_dist = cur_dist[order(apply(cur_dist,MARGIN = 1,FUN = mean),decreasing = T),]
      singles  = rep(TRUE,length(cur_dist[,1]))
      couples  = repmat(NaN,length(cur_dist[,1]),3)
      for (kk in 1:length(cur_dist[,1])) {
        couples[kk,1]    = kk
        cur_ranks        = rank(cur_dist[kk,])
        cur_av_rnks      = cur_ranks[singles]
        min_rank         = min(cur_av_rnks)
        partner          = which(cur_ranks == min_rank)[1]
        couples[kk,2]    = partner
        couples[kk,3]    = cur_dist[kk,partner]
        singles[partner] = FALSE
      }
      
      # which one is the most distanced couple?
      ct = ct + 1
      couple_to_kill = which(max(couples[,3]) == couples[,3])
      HC_kill[ct]    = rownames(cur_dist)[couples[couple_to_kill,1]]
      PG_kill[ct]    = colnames(cur_dist)[couples[couple_to_kill,2]]
      
      # actually kill the couple
      # in case of size=1 for cur_names_dom_narrowed, we need to recover names and rownames
      names_HC                 = names(cur_df_HC_vars)
      names_PG                 = names(cur_df_PG_vars)
      rnames_HC                = rownames(cur_df_HC_vars)
      rnames_PG                = rownames(cur_df_PG_vars)
      cur_df_HC_vars           = as.data.frame(cur_df_HC_vars[-which(HC_kill[ct] == rownames(cur_df_HC_vars)),])
      cur_df_PG_vars           = as.data.frame(cur_df_PG_vars[-which(PG_kill[ct] == rownames(cur_df_PG_vars)),])
      rnames_HC                = rnames_HC[-which(HC_kill[ct] == rnames_HC)]
      rnames_PG                = rnames_PG[-which(PG_kill[ct] == rnames_PG)]
      names(cur_df_HC_vars)    = names_HC
      names(cur_df_PG_vars)    = names_PG
      rownames(cur_df_HC_vars) = rnames_HC
      rownames(cur_df_PG_vars) = rnames_PG
      
      # test improvement in matching
      cur_df     = cur_df[!row.names(cur_df) %in% c(HC_kill,PG_kill),]
      dfs[[ii]]  = cur_df
      res        = agk.perform.matching.tests(dfs,cur_groups,cur_matching,path_mtc,0,cur_names)
      score      = res[[ii]]$ok_or_checkneeded[res[[ii]]$cur_var %in% cur_names_dom_narrowed]
      contn_elim = any(score == 'NEEDS_CHECK')
    }
    # saving the dropped PGs and HCs
    cur_drop_list = list()
    cur_drop_list$dropped_HC_matching = HC_kill
    cur_drop_list$dropped_PG_matching = PG_kill
    drop_lists[[ii]]                  = cur_drop_list
    
  }
  return(list(dfs = dfs, dropped_HCs_PGs = drop_lists))
}