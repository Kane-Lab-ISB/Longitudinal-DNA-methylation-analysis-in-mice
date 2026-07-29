#colors
#female  "#AA4466" male "#4477AA"
#time points BL-T5 BL = "#332288", T2 = "#AAAA00", T3 = "#44AA99", T4 = "#882255", T5 = "#88CCEE"
#manhattan plot  "#114477" "#44AA99"  
#four categories cpg context; gene context; XCI status, 
#"#DDCC77", "#CC6677", "#AA4499", "#882255"  
#three categories aicpg/dmr  accpg/dmr dualcpg/dmr
#frail "#762A83" age"#44AA99" dual "#DDAA33" other "#999933"
#visualization background "#f5f0f0"

#dmr_fic "#762A83" 
#dmr_delta "#117777" 
#delta_fic "#AA7744" 
#delta_fif "#DDAA33"
#delta_delta "#88CC88"

#"#AA4466", "#4477AA", "#8844AA", "#117777", "#66AA55", "#AA7744", "#DDAA33"

get_fit2_res = function(metadata, design, values) {
  dat_limma = values %>% filter(rownames(.) %in% metadata$id)
  counts_excl = dat_limma %>% t() %>% as.data.frame() 
  corfit = duplicateCorrelation(counts_excl, design, block = metadata$Mouse.ID)
  fit_ns = lmFit(counts_excl, design, block = metadata$Mouse.ID, correlation = corfit$consensus)
  fit2_ns = eBayes(fit_ns, trend = TRUE, robust = TRUE)
  return(fit2_ns)
} ## change this in the end

get_dmr_modify_DMRcate = function(cpg_fc) {
  C = 2
  lambda = 1000
  consec = FALSE
  lag = lambda
  if (is.null(cpg_fc$`F`)) {
    fc_change_add = cpg_fc %>%
      mutate(is.sig = ifelse(adj.P.Val < 0.05, TRUE, FALSE)) %>% 
      mutate(weights = abs(t)) 
  } else {
    fc_change_add = cpg_fc %>%
      mutate(is.sig = ifelse(adj.P.Val < 0.05, TRUE, FALSE)) %>% 
      mutate(weights = sqrt(`F`)) 
  }
  object = data.frame(ID = fc_change_add$IlmnID,
                      weights = fc_change_add$weights,
                      CHR = as.character(fc_change_add$chromo),
                      pos = as.numeric(fc_change_add$MAPINFO),
                      indfdr = fc_change_add$adj.P.Val,
                      diff = rep(0, length(fc_change_add$IlmnID)),
                      is.sig = fc_change_add$is.sig)
  object$CHR = factor(object$CHR, levels = chrOrder_DNAm)
  object = object[order(object$CHR, object$pos), ]
  fitted = lapply(chrOrder_DNAm, function(x) {
    chromosome = object[object$CHR == x, ]
    pos = chromosome$pos
    sigma = lambda/C
    lag = lambda
    beta = chromosome$weights
    df = 1
    X2 = beta^2
    pvalue = KernelTest(pos = pos, X2 = X2, lambda = sigma, df = df)
    chromosome$raw = pvalue
    return(chromosome)  
  })
  object = rbind.fill(fitted)
  object$fdr = p.adjust(object$raw, method = "BH")
  nsig = sum(object$is.sig)
  pcutoff = sort(object$fdr)[nsig]
  object$sig = (object$fdr <= pcutoff)
  K = cpg_fc %>% filter(adj.P.Val < 0.05) %>% nrow()
  k = K - 1
  chr.N = as.character(object$CHR)
  pos.N = object$pos
  sig.N = object$sig
  N = length(sig.N)
  n.K = which(sig.N)
  K = length(n.K)
  pos.K = pos.N[n.K]
  chr.K = chr.N[n.K]
  jump_chr.k = (chr.K[-1] != chr.K[-K])
  jump_pos.k = (diff(pos.K) > lag)
  jump.k = (jump_chr.k | jump_pos.k)
  ksegments.A2 = Segment(jump.k)
  A = nrow(ksegments.A2)
  kstart.A = ksegments.A2[, "start"]
  kend.A = ksegments.A2[, "end"]
  realpos.K = pos.K
  start.A = realpos.K[kstart.A]
  end.A = realpos.K[kend.A]
  chr.A = chr.K[kstart.A]
  fmt = "%s:%1d-%1d"
  coord.A = sprintf(fmt, chr.A, start.A, end.A)
  nstart.A = n.K[kstart.A]
  nend.A = n.K[kend.A]
  width.A = nend.A + 1 - nstart.A
  a.Z = rep(seq(A), width.A)
  fn = function(a) seq(from = nstart.A[a], to = nend.A[a])
  l.listA = lapply(seq(A), fn)
  n.Z = unlist(l.listA)
  region.N = rep(NA_integer_, N)
  region.N[n.Z] = a.Z
  levels = seq(A)
  region.N = factor(region.N, levels = levels)
  no_cpg.A = c(table(region.N))
  REGIONSTAT = function(field, fn) {
    x.N = object[[field]]
    x.R = tapply(x.N, region.N, fn)
    c(x.R)
  }
  fn_Stouffer = function(x) pnorm(sum(qnorm(x))/sqrt(length(x)))
  fn_HMFDR = function (x) 1/mean(1/x)
  fn_Fisher = function (x) pchisq((sum(log(x))*-2), df=length(x)*2, lower.tail=FALSE)
  fn_max = function(x) x[which.max(abs(x))]
  results = data.frame(
    coord = coord.A,
    no.cpgs = no_cpg.A,
    chr = chr.A,
    start = start.A,
    end = end.A,
    min_smoothed_fdr = REGIONSTAT("fdr", min),
    Stouffer = REGIONSTAT("indfdr", fn_Stouffer),
    HMFDR = REGIONSTAT("indfdr", fn_HMFDR),
    Fisher = REGIONSTAT("indfdr", fn_Fisher),
    maxdiff = REGIONSTAT("diff", fn_max),
    meandiff = REGIONSTAT("diff", mean),
    row.names = seq(A),
    stringsAsFactors = FALSE)
  return(results)
}

select_dmp_for_dmr = function(dmr_lst) {
  all_dmp = lapply(1:nrow(dmr_lst), function(x) {
    chr_dmr = dmr_lst[x, ]$chr
    selected_dmp = cpg_loci %>%
      filter(IlmnID %in% colnames(methy_dat_dis_M)) %>%
      filter(chromo == chr_dmr) %>%
      filter(MAPINFO >= dmr_lst[x, ]$start & MAPINFO <= dmr_lst[x, ]$end) %>%
      mutate(dmr = dmr_lst[x, ]$coord) %>%
      dplyr::select(IlmnID, Name, chromo, MAPINFO, dmr)
    return(selected_dmp)
  })
  dmp_w_dmr = Reduce(rbind, all_dmp) %>%
    as.data.frame()
  return(dmp_w_dmr)
}

prob_enriched_genes = function(cpg_list) {
  res = testEnrichment(cpg_list, 
                       buildGeneDBs(cpg_list, max_distance = 10000, platform = "MM285"),
                       platform = "MM285") %>%
    filter(FDR < 0.05) %>%
    filter(as.numeric(overlap)/as.numeric(nD) > 0.05)
  return(res)
}

prob_enriched_path = function(cpg_list) {
  res = testGO(cpg_list, platform = "MM285", organism = "mmusculus")
  result = res$result %>% 
    filter(source == "KEGG") %>% 
    filter(significant == "TRUE") %>% 
    mutate(prop = intersection_size/term_size) %>% arrange(desc(p_value))
  return(result)
}

test_rdrop = function(loci_dmr, expr_tbl) {
  rdrop_tbl = lapply(unique(loci_dmr$dmr), function(x) {
    single_dmr = loci_dmr %>%
      filter(dmr == x)
    test_matrix = expr_tbl %>%
      dplyr::select(single_dmr$IlmnID)
    rdop_num = sapply(1:ncol(test_matrix), function(y) {
      data_y = test_matrix[, y]
      data_no_y = test_matrix[, -y] %>% rowMeans()
      num = cor(data_y, data_no_y, method = "spearman")
    })
    single_dmr_rdrop = single_dmr %>%
      mutate(rdrop = rdop_num)
    return(single_dmr_rdrop)
  }) %>%
    Reduce(rbind, .) %>%
    as.data.frame()
  return(rdrop_tbl)  
} 

prepostfuture_data_merge = function(tp_group, assoc_df, fea_lst) { #function to get time point pair
  merge_dat = lapply(1:nrow(tp_group), function(x) {
    dat = get_prepostfuture_dat(tp_group[[x, 1]], tp_group[[x, 2]], tp_group[[x, 3]],
                          assoc_df, fea_lst) %>%
      mutate(level = paste0(tp_group[[x, 3]], "_", tp_group[[x, 2]], "_", tp_group[[x, 1]]))
    return(dat)
  }) %>%
    Reduce(rbind, .) %>%
    as.data.frame()
  return(merge_dat)
}

prepost_data_merge = function(tp_group, assoc_df, fea_lst) { #function to get time point pair
  merge_dat = lapply(1:nrow(tp_group), function(x) {
    dat = get_prepost_dat(tp_group[[x, 1]], tp_group[[x, 2]], 
                          assoc_df, fea_lst) %>%
      mutate(level = paste0(tp_group[[x, 2]], "_", tp_group[[x, 1]]))
    return(dat)
  }) %>%
    Reduce(rbind, .) %>%
    as.data.frame()
  return(merge_dat)
}

get_prepost_dat = function(pre, post, assoc_df, fea_lst) {#function to get time point pair
  pre_data = assoc_df %>% filter(Time.label == pre)
  post_data = assoc_df %>% filter(Time.label == post)
  prepost_data = merge(pre_data, post_data, by = "mouse_id") %>%
    as.data.frame() %>%
    mutate(age_change = `age.y` - `age.x`,
           fi_change = `logFI.y`/`logFI.x`)
  deltaMeth = lapply(fea_lst, function(x){
    change = prepost_data[, paste0(x, ".y")] - prepost_data[, paste0(x, ".x")]
    return(change)
  }) 
  deltaMeth_df = Reduce(cbind, deltaMeth) %>%
    as.data.frame() %>%
    `colnames<-`(paste0("delta_", fea_lst))
  prepost_data_out = cbind(prepost_data, deltaMeth_df) %>%
    as.data.frame()
}

get_prepostfuture_dat = function(pre, post, future, assoc_df, fea_lst) {#function to get time point pair
  pre_data = assoc_df %>% filter(Time.label == pre)
  post_data = assoc_df %>% filter(Time.label == post)
  future_data = assoc_df %>% filter(Time.label == future)
  prepost_data = merge(pre_data, post_data, by = "mouse_id") %>%
    merge(., future_data, by = "mouse_id") %>%
    as.data.frame() %>%
    mutate(age_change1 = `age.y` - `age.x`,
           age_change2 = `age` - `age.y`,
           fi_change = `logFI`/`logFI.y`)
  deltaMeth = lapply(fea_lst, function(x){
    change = prepost_data[, paste0(x, ".y")] - prepost_data[, paste0(x, ".x")]
    return(change)
  }) 
  deltaMeth_df = Reduce(cbind, deltaMeth) %>%
    as.data.frame() %>%
    `colnames<-`(paste0("delta_", fea_lst))
  prepost_data_out = cbind(prepost_data, deltaMeth_df) %>%
    as.data.frame()
}

dmr_assoc_frailty_outcome = function(frail_assoc_df, frail_prepost_assoc_df, frail_prepostfuture_assoc_df) {
  dmr_assoc = lapply(1:925, function(x) {
    loci_lmer = fidmrs_cpg_combined %>% filter(dmr == unique(fidmrs_cpg_combined$dmr)[[x]])
    loci_names1 = loci_lmer$IlmnID
    loci_names2 = paste0(loci_lmer$IlmnID, ".x")
    loci_names3 = paste0("delta_", loci_lmer$IlmnID)
    expr_merged1 = frail_assoc_df %>%
      dplyr::select(all_of(loci_names1), id, mouse_id, age, logFI, sex) %>%
      melt(., id = c("id", "mouse_id", "age", "logFI", "sex")) %>%
      as.data.frame()
    res1 = lmer(logFI ~ value + age + as.factor(sex) + (1 | mouse_id) + (1 | variable), 
                data = expr_merged1)
    res2 = lmer(logFI ~ value*age + as.factor(sex) + (1 | mouse_id) + (1 | variable), 
                data = expr_merged1)
    expr_merged2 = frail_prepost_assoc_df %>%
      dplyr::select(all_of(loci_names2), id.x, mouse_id, age.x, age_change, logFI.y, fi_change, sex.x) %>%
      melt(., id = c("id.x", "mouse_id", "age.x", "age_change", "logFI.y", "fi_change", "sex.x")) %>%
      as.data.frame()
    res3 = lmer(logFI.y ~ value + age.x + as.factor(sex.x) + age_change + (1 | mouse_id) + (1 | variable), 
                data = expr_merged2)
    res4 = lmer(logFI.y ~ value*age.x + as.factor(sex.x) +  age_change + (1 | mouse_id) + (1 | variable), 
                data = expr_merged2)
    res5 = lmer(fi_change ~ value + age.x + as.factor(sex.x) + age_change + (1 | mouse_id) + (1 | variable), 
                data = expr_merged2)
    res6 = lmer(fi_change ~ value*age.x + as.factor(sex.x) + age_change + (1 | mouse_id) + (1 | variable), 
                data = expr_merged2)
    expr_merged3 = frail_prepost_assoc_df %>%
      dplyr::select(all_of(loci_names3), id.x, mouse_id, age.y, age_change, logFI.y, sex.x) %>%
      melt(., id = c("id.x", "mouse_id", "age.y", "age_change", "logFI.y", "sex.x")) %>%
      as.data.frame()
    res7 = lmer(logFI.y ~ value + age.y + as.factor(sex.x) + age_change + (1 | mouse_id) + (1 | variable), 
                data = expr_merged3)
    res8 = lmer(logFI.y ~ value*age.y + as.factor(sex.x) +  age_change + (1 | mouse_id) + (1 | variable), 
                data = expr_merged3)
    expr_merged4 = frail_prepostfuture_assoc_df %>%
      dplyr::select(all_of(loci_names3), id.x, mouse_id, age.y, age_change1, age_change2, logFI, fi_change, sex) %>%
      melt(., id = c("id.x", "mouse_id", "age.y", "age_change1", "age_change2", "logFI", "fi_change", "sex")) %>%
      as.data.frame()
    res9 = lmer(logFI ~ value + age.y + as.factor(sex) + age_change1 + age_change2 + (1 | mouse_id) + (1 | variable),
                data = expr_merged4)
    res10 = lmer(logFI ~ value*age.y + as.factor(sex) +  age_change1 + age_change2 + (1 | mouse_id) + (1 | variable), 
                 data = expr_merged4)
    res11 = lmer(fi_change ~ value + age.y + as.factor(sex) + age_change1 + age_change2 + (1 | mouse_id) + (1 | variable), 
                 data = expr_merged4)
    res12 = lmer(fi_change ~ value*age.y + as.factor(sex) + age_change1 + age_change2 + (1 | mouse_id) + (1 | variable), 
                 data = expr_merged4)
    return(c(res1, res2, res3, res4, res5, res6, res7, res8, res9, res10, res11, res12))
  })
}

dmr_assoc_frailty_outcome_sex = function(frail_assoc_df, frail_prepost_assoc_df, frail_prepostfuture_assoc_df) {
  dmr_assoc = lapply(1:925, function(x) {
    loci_lmer = fidmrs_cpg_combined %>% filter(dmr == unique(fidmrs_cpg_combined$dmr)[[x]])
    loci_names1 = loci_lmer$IlmnID
    loci_names2 = paste0(loci_lmer$IlmnID, ".x")
    loci_names3 = paste0("delta_", loci_lmer$IlmnID)
    expr_merged1 = frail_assoc_df %>%
      dplyr::select(all_of(loci_names1), id, mouse_id, age, logFI) %>%
      melt(., id = c("id", "mouse_id", "age", "logFI")) %>%
      as.data.frame()
    res1 = lmer(logFI ~ value + age + (1 | mouse_id) + (1 | variable), 
                data = expr_merged1)
    res2 = lmer(logFI ~ value*age + (1 | mouse_id) + (1 | variable), 
                data = expr_merged1)
    expr_merged2 = frail_prepost_assoc_df %>%
      dplyr::select(all_of(loci_names2), id.x, mouse_id, age.x, age_change, logFI.y, fi_change) %>%
      melt(., id = c("id.x", "mouse_id", "age.x", "age_change", "logFI.y", "fi_change")) %>%
      as.data.frame()
    res3 = lmer(logFI.y ~ value + age.x + age_change + (1 | mouse_id) + (1 | variable), 
                data = expr_merged2)
    res4 = lmer(logFI.y ~ value*age.x + age_change + (1 | mouse_id) + (1 | variable), 
                data = expr_merged2)
    res5 = lmer(fi_change ~ value + age.x + age_change + (1 | mouse_id) + (1 | variable), 
                data = expr_merged2)
    res6 = lmer(fi_change ~ value*age.x + age_change + (1 | mouse_id) + (1 | variable), 
                data = expr_merged2)
    expr_merged3 = frail_prepost_assoc_df %>%
      dplyr::select(all_of(loci_names3), id.x, mouse_id, age.y, age_change, logFI.y) %>%
      melt(., id = c("id.x", "mouse_id", "age.y", "age_change", "logFI.y")) %>%
      as.data.frame()
    res7 = lmer(logFI.y ~ value + age.y + age_change + (1 | mouse_id) + (1 | variable), 
                data = expr_merged3)
    res8 = lmer(logFI.y ~ value*age.y + age_change + (1 | mouse_id) + (1 | variable), 
                data = expr_merged3)
    expr_merged4 = frail_prepostfuture_assoc_df %>%
      dplyr::select(all_of(loci_names3), id.x, mouse_id, age.y, age_change1, age_change2, logFI, fi_change) %>%
      melt(., id = c("id.x", "mouse_id", "age.y", "age_change1", "age_change2", "logFI", "fi_change")) %>%
      as.data.frame()
    res9 = lmer(logFI ~ value + age.y + age_change1 + age_change2 + (1 | mouse_id) + (1 | variable),
                data = expr_merged4)
    res10 = lmer(logFI ~ value*age.y + age_change1 + age_change2 + (1 | mouse_id) + (1 | variable), 
                 data = expr_merged4)
    res11 = lmer(fi_change ~ value + age.y + age_change1 + age_change2 + (1 | mouse_id) + (1 | variable), 
                 data = expr_merged4)
    res12 = lmer(fi_change ~ value*age.y + age_change1 + age_change2 + (1 | mouse_id) + (1 | variable), 
                 data = expr_merged4)
    return(c(res1, res2, res3, res4, res5, res6, res7, res8, res9, res10, res11, res12))
  })
}
