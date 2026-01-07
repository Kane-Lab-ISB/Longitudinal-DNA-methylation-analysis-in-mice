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
}

tt_aging_lst = function(fit2) {
  cont_mix = makeContrasts((sexfemale.X3 + sexmale.X3)/2 - (sexfemale.X2 + sexmale.X2)/2,
                           (sexfemale.X2 + sexmale.X2)/2 - (sexfemale.X1 + sexmale.X1)/2,
                           levels = colnames(coef(fit2)))
  cont_female = makeContrasts(sexfemale.X2 - sexfemale.X1,
                              sexfemale.X3 - sexfemale.X2,
                              levels = colnames(coef(fit2)))
  cont_male = makeContrasts(sexmale.X2 - sexmale.X1,
                            sexmale.X3 - sexmale.X2,
                            levels = colnames(coef(fit2)))
  cont_sex = makeContrasts(sexfemale.X1 - sexmale.X1,
                           sexfemale.X2 - sexmale.X2,
                           sexfemale.X3 - sexmale.X3,
                           levels = colnames(coef(fit2)))
  cont_lst = list(cont_mix, cont_female, cont_male, cont_sex)
  tt_global_lst = lapply(cont_lst, function(y) {
    fit_cont = contrasts.fit(fit2, y)
    fit_cont = eBayes(fit_cont)
    tt_cont = topTable(fit_cont, n = Inf, adjust.method = "BH") 
    return(tt_cont)
  })
  return(tt_global_lst)
}

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

methy_prob_loci = read.csv(file = "~/Desktop/DNA methylation/MouseMethylation-12v1-0_A2.csv")

get_DMRs = function(values, tt) {
  cpg_loci = methy_prob_loci %>%
    filter(IlmnID %in% colnames(values)) %>%
    mutate(chromo = paste0("chr", CHR))
  cpg_fc = merge(tt %>% mutate(IlmnID = rownames(.)), cpg_loci, by = "IlmnID")
  DMR_res = get_dmr_modify_DMRcate(cpg_fc) %>%
    filter(no.cpgs >= 5 & Stouffer < 0.05)
  return(DMR_res)
}

manhattan_plot_ewas = function(tt_tbl) {
  cpg_pos = merge(tt_tbl %>% mutate(IlmnID = rownames(.)),
                  cpg_loci,
                  by = "IlmnID") %>% 
    as.data.frame() %>%
    dplyr::select(IlmnID, chromo, MAPINFO, P.Value, adj.P.Val)
  cpg_pos_cumpos = merge(cpg_pos, chromo_cumpos, by = "chromo") %>%
    mutate(bpcum = MAPINFO + position_x) %>%
    dplyr::select(chromo, IlmnID, P.Value, adj.P.Val, MAPINFO, position_x, bpcum)
  cpg_pos_cumpos$chromo = factor(cpg_pos_cumpos$chromo, levels = chrOrder_DNAm)
  nondmp_pos_cumpos = cpg_pos_cumpos %>% filter(adj.P.Val >= 0.05)
  dmp_pos_cumpos = cpg_pos_cumpos %>% filter(adj.P.Val < 0.05)
  x_axis_ticks = cpg_pos_cumpos %>%
    dplyr::group_by(chromo) %>%
    dplyr::summarise(center = (max(bpcum) + min(bpcum))/2)
  plot = ggplot() +
    geom_point(data = nondmp_pos_cumpos,
               aes(x = bpcum,
                   y = -log10(P.Value),
                   color = as.factor(chromo)), size = 0.2) + 
    scale_color_manual(values = rep(c("grey65", "grey80"), 21)) +
    guides(color = FALSE) + 
    new_scale_color() + 
    geom_point(data = dmp_pos_cumpos,
               aes(x = bpcum,
                   y = -log10(P.Value), 
                   color = as.factor(chromo)), size = 0.4) + 
    scale_color_manual(values = rep(c("#004488", "#44AA99"), 22)) + 
    scale_x_continuous(label = x_axis_ticks$chromo, breaks = x_axis_ticks$center) +
    labs(x = "", y = "-log10(p-value)", color = "") +
    guides(color = FALSE) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 12, face = "bold"))
  return(plot)
}

manhattan_plot_dmr = function(tt_tbl, dmrs_tbl, loci_dmrs_tbl) {
  cpg_pos = merge(tt_tbl %>% mutate(IlmnID = rownames(.)),
                  cpg_loci,
                  by = "IlmnID") %>% 
    as.data.frame() %>%
    dplyr::select(IlmnID, chromo, MAPINFO, P.Value, adj.P.Val)
  cpg_pos_cumpos = merge(cpg_pos, chromo_cumpos, by = "chromo") %>%
    mutate(bpcum = MAPINFO + position_x) %>%
    dplyr::select(chromo, IlmnID, P.Value, adj.P.Val, MAPINFO, position_x, bpcum)
  cpg_pos_cumpos$chromo = factor(cpg_pos_cumpos$chromo, levels = chrOrder_DNAm)
  dmr_pos = dmrs_tbl %>%
    mutate(MAPINFO = (start + end)/2,
           IlmnID = coord,
           chromo = chr,
           adj.P.Val = Stouffer) %>%
    dplyr::select(IlmnID, chromo, MAPINFO, adj.P.Val) %>%
    mutate(type = rep("dmr", nrow(.)))
  dmr_pos_cumpos = merge(dmr_pos, chromo_cumpos, by = "chromo") %>%
    mutate(bpcum = MAPINFO + position_x)
  nondmp_pos_cumpos = cpg_pos_cumpos %>% filter(adj.P.Val >= 0.05) %>% filter(IlmnID %in% loci_dmrs_tbl$IlmnID)
  dmp_pos_cumpos = cpg_pos_cumpos %>% filter(adj.P.Val < 0.05) %>% filter(IlmnID %in% loci_dmrs_tbl$IlmnID)
  x_axis_ticks = cpg_pos_cumpos %>%
    dplyr::group_by(chromo) %>%
    dplyr::summarise(center = (max(bpcum) + min(bpcum))/2)
  plot = ggplot() +
    geom_point(data = nondmp_pos_cumpos,
               aes(x = bpcum,
                   y = -log10(P.Value),
                   color = as.factor(chromo)), size = 0.2) + 
    scale_color_manual(values = rep(c("grey65", "grey80"), 21)) +
    guides(color = FALSE) + 
    new_scale_color() + 
    geom_point(data = dmp_pos_cumpos,
               aes(x = bpcum,
                   y = -log10(P.Value), 
                   color = as.factor(chromo)), size = 0.4) + 
    scale_color_manual(values = rep(c("#BADE28FF", "#5BC863FF"), 22)) + 
    geom_point(data = dmr_pos_cumpos,
               aes(x = bpcum,
                   y = -log10(adj.P.Val)),
               color = "#443983FF", shape = 2) +
    scale_x_continuous(label = x_axis_ticks$chromo, breaks = x_axis_ticks$center) +
    labs(x = "", y = "-log10(p-value)", color = "") +
    guides(color = FALSE) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 12, face = "bold"))
  return(plot)
}


chrOrder_DNAm = c(paste0("chr", 1:19), "chrX", "chrY", "chrMT")
cpg_loci$chromo = factor(cpg_loci$chromo, levels = chrOrder_DNAm)
chromo_cumpos = cpg_loci %>%
  group_by(chromo) %>%
  dplyr::summarise(chr_len = max(MAPINFO)) %>%
  mutate(chr_len = as.numeric(chr_len)) %>%
  mutate(position_x = cumsum(chr_len) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  as.data.frame()

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

enriched_gene_plot = function(gene_output) {
  top25 = gene_output %>%  
    dplyr::slice(1:20) %>% 
    mutate(ratio = overlap/nD) %>% 
    arrange(desc(p.value))
  top25$gene_name = factor(top25$gene_name, levels = top25$gene_name)
  plot = ggplot(data = top25,
                aes(x = -log10.p.value,
                    y = gene_name,
                    fill = ratio)) +
    geom_bar(stat = "identity", width = 0.75) +
    scale_fill_viridis() +
    labs(y = "") +
    theme(axis.text.y = element_text(size = 9, face = "bold.italic"))
  return(plot)
}

kegg_path_plot = function(path_output) {
  top20 = path_output %>% arrange(p_value) %>% 
    dplyr::slice(1:20) %>% arrange(desc(p_value))
  top20$term_name = factor(top20$term_name, levels = top20$term_name)
  plot = ggplot(data = top20,
                aes(x = intersection_size,
                    y = term_name,
                    color = -log10(p_value),
                    size = prop)) +
    geom_point()+
    scale_color_viridis(direction = -1) +
    theme(axis.text.y = element_text(size = 8, face = "bold")) +
    labs(y = "")
  return(plot)
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

get_coeff = function(fea_lst, lmer_lst) { # get the coefficient table for association study res
  coeff_tbl_merge = sapply(1:length(fea_lst), function(x) {
    coeff_tbl = summary(lmer_lst[[x]])$coefficients
    return(c(coeff_tbl[2, 1], coeff_tbl[2, 2], coeff_tbl[2, 5]))
  }) %>% t() %>% 
    as.data.frame() %>%
    mutate(cpg = unlist(fea_lst)) %>%
    `colnames<-`(c("coeff", "se", "p", "cpg")) %>%
    mutate(adjp = p.adjust(p, method = "fdr")) %>% 
    filter(adjp < 0.05)
  return(coeff_tbl_merge)
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

get_freq_dmr = function(dmr_fic_df, dmr_fif_df, dmr_deltafi_df,
                        deltadmr_fic_df, deltadmr_fif_df, deltadmr_deltafi_df) {
  dmr_assoc_list = Reduce(union, list(dmr_fic_df$dmr,
                                      dmr_fif_df$dmr,
                                      dmr_deltafi_df$dmr,
                                      deltadmr_fic_df$dmr,
                                      deltadmr_fif_df$dmr,
                                      deltadmr_deltafi_df$dmr))
  dmr_assoc_freq = data.frame(dmr_fic = ifelse(dmr_assoc_list %in% dmr_fic_df$dmr, 1, 0),
                              dmr_fif = ifelse(dmr_assoc_list %in% dmr_fif_df$dmr, 1, 0),
                              dmr_deltafi = ifelse(dmr_assoc_list %in% dmr_deltafi_df$dmr, 1, 0),
                              deltadmr_fic = ifelse(dmr_assoc_list %in% deltadmr_fic_df$dmr, 1, 0),
                              deltadmr_fif = ifelse(dmr_assoc_list %in% deltadmr_fif_df$dmr, 1, 0),
                              deltadmr_deltafi = ifelse(dmr_assoc_list %in% deltadmr_deltafi_df$dmr, 1, 0)) %>%
    mutate(num = rowSums(.),
           dmr = dmr_assoc_list) %>%
    filter(num >= 2) 
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