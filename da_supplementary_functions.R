# ex. there are two 'f__', both from different orders.
update_taxonomy = function(ps){
  tax_info = ps@tax_table %>% as.matrix() %>% as.data.frame() %>% 
    rownames_to_column('OTU')
  tax_info = tax_info[,which(!is.na(tax_info[1,]))]
  
  info_column = ncol(tax_info)-1 # Initialize the column used to add tax info
  while(length(unique(tax_info[,ncol(tax_info)]))<length(tax_info[,ncol(tax_info)]) & info_column>0){
    temp = tax_info[,ncol(tax_info)]
    temp_dup = tax_info[duplicated(temp),] # Select only rows that have duplicates in final col
    temp_dup[,ncol(temp_dup)] = paste(temp_dup[,info_column],temp_dup[,ncol(temp_dup)],sep='.')
    tax_info[duplicated(temp),] = temp_dup
    info_column = info_column-1
  }
  if(length(unique(tax_info[,ncol(tax_info)]))<length(tax_info[,ncol(tax_info)])){
    print('There are duplicate taxonomic identities that could not be resolved using upper level taxonomy. Please inspect the output taxonomic table.')
    print(length(unique(tax_info[,ncol(tax_info)])))
    print(length(tax_info[,ncol(tax_info)]))
    return(tax_info)
    }
  
  return(tax_info)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tidy_ancombc = function(phyloseq, formula, p_adj_method = "fdr", prv_cut = 0.1, abun_cut = 0.0001,
                        lib_cut = 0, group = NULL, struc_zero = FALSE, neg_lb = FALSE, tax_level = NULL,
                        tol = 1e-05, max_iter = 100, conserve = FALSE, alpha = 0.05,
                        global = FALSE,n_cl = 1){
  
  # phyloseq=x;formula='TPN_control';p_adj_method = "fdr"; prv_cut = 0.1;abun_cut = 0;
  # lib_cut = 0; group = 'TPN_control'; struc_zero = T; neg_lb = F;
  # tol = 1e-05; max_iter = 100; conserve = TRUE; alpha = 0.05;
  # global = FALSE;n_cl=1;tax_level='Species'
  
  # Set seed and run ancom-bc
  set.seed(421)
  ancom.sp = ancombc(data=phyloseq, formula = formula, p_adj_method=p_adj_method, prv_cut=prv_cut,
                     lib_cut=lib_cut, group=group, struc_zero=struc_zero, neg_lb = neg_lb,tax_level=tax_level,
                     tol=tol, max_iter=max_iter, conserve=conserve, alpha=alpha,global=global,n_cl=n_cl)
  print('ANCOM done')
  # Determine key parameters for other functions based on input data
  ## Taxonomic level
  t = which(!is.na(phyloseq@tax_table[1,])) %>% max # select lowest level
  t = colnames(phyloseq@tax_table)[t] # Name of tax level
  t2 = str_sub(t,1,1) %>% str_to_lower() # ex. 's' for Species
  
  # Format ancom-bc results
  print('Starting formatting')
  if(is.null(tax_level)) tlevel=colnames(ps@tax_table)[1] else tlevel = tax_level
  ancom.res = format_ancom_results(ancom.sp,phyloseq %>% tax_glom(tlevel),abun_cut,tlevel)
  
  if(is.null(group)){group=NA}
  
  # Settings used for the analysis
  settings = data.frame(Parameter = c('tax_level', 'tax_prefix', 'formula', 'p_adj_method', 'prv_cut',
                                      'abun_cut',
                                      'lib_cut', 'group', 'struc_zero', 'neg_lb', 'tol', 'max_iter','conserve', 'alpha', 'global'),
                        Value = c(t, t2, formula, p_adj_method, prv_cut, abun_cut,
                                  lib_cut, group, struc_zero, neg_lb, tol, max_iter,conserve, alpha, global))
  
  return(list(results = ancom.res,input_ps = phyloseq,settings = settings))
}
#~~~~~~~~~~

tidy_aldex2 =  function(phyloseq,formula,aldex_type = 'kw',prv_cut = 0.10,abun_cut=0.0001,
                       mc.samples=128,denom="all",useMC=FALSE,decimals = F){
  
  # phyloseq = ps %>% tax_glom('Phylum'); formula = 'Status';
  # aldex_type = 'kw'; prv_cut = 1/3;abun_cut=1e-4;mc.samples=128;denom='all';useMC=F;decimals=F

  ## Taxonomic level
  t = which(!is.na(phyloseq@tax_table[1,])) %>% max # select lowest level
  t = colnames(phyloseq@tax_table)[t] # Name of tax level
  t2 = str_sub(t,1,1) %>% str_to_lower() # ex. 's' for Species
  
  # Extract variables
  f = formula %>% str_split('[+]') %>% lapply(str_trim) %>% unlist()
  
  s = phyloseq@sam_data %>% as.matrix() %>% as.data.frame()
  # Convert to numeric where applicable
  # Note: if the first value in a column is NA, then it won't convert properly.
  w = which(!is.na(as.numeric(s[1,]))) %>% suppressWarnings()
  s[, w] <- lapply(s[, w, drop = FALSE], function(x) as.numeric(x))
  o = phyloseq@otu_table %>% as.matrix() %>% as.data.frame()
  o = o %>% dplyr::select(all_of(rownames(s))) 
  
  if(aldex_type=='kw'){ 
    m = s[[f]] 
  } else if(aldex_type == 'glm'){ 
    m=eval(parse(text=paste("model.matrix(~",formula,",data=s)",sep='')))
  } else if(aldex_type == 'spearman'){ 
    m=eval(parse(text=paste("model.matrix(~",formula,",data=s)",sep='')))
  }
  if(!exists('m')){ print('Error: aldex_type must be kw, spearman, or glm.'); break }
  
  # if 'kw', need to extract the second level of the variable for column names
  if(aldex_type=='kw'){ 
    v = factor(m);v2 = length(levels(v))
    if(v2>2){ 
      res.var=f 
      print('Variable has more than two levels; colnames will use variable as prefix without level information.')
    } else {
      res.var=paste(f,levels(v)[2],sep='')
    }
  }
  
  # Run ALDEX2
  set.seed(421)
  if(decimals==T){
    o2 = o*(1/min(o[o!=0])) # Smallest nonzero value transformed to 1 read (so that rounding doesn't remove significant variance)
    x = aldex.clr(round(o2,0),m,mc.samples = mc.samples, denom = denom, verbose=F)
  } else {
    x = aldex.clr(o,m,mc.samples = mc.samples, denom = denom)
  }
  if(aldex_type=='kw'){ 
    aldex.sp = aldex.kw(x) 
    x.effect <- aldex.effect(x, CI=T, verbose=FALSE, paired.test=FALSE)
    aldex.sp = merge(aldex.sp,x.effect,by='row.names',all=T)
    aldex.res = format_aldex_results(aldex.sp,phyloseq,x,res.var,min_prevalence = prv_cut,
                                     min_abundance=abun_cut, aldex_type=aldex_type)
  } else if(aldex_type=='glm'){ 
    aldex.sp = aldex.glm(x) 
    x.effect <- aldex.glm.effect(x, CI=T, verbose=FALSE)
    aldex.sp = merge(aldex.sp,x.effect,by='row.names',all=T)
    aldex.res = format_aldex_results(aldex.sp,phyloseq,x,min_prevalence = prv_cut,
                                     min_abundance = abun_cut, aldex_type=aldex_type)
  } else if(aldex_type=='spearman'){
    aldex.sp = aldex.corr(x,s$mindscore) %>% suppressWarnings()
    aldex.res = format_aldex_results(aldex.sp,phyloseq,x,res.var=f,min_prevalence = prv_cut,
                                     min_abundance = abun_cut, aldex_type=aldex_type)
  }
  
  settings = data.frame(Parameter = c('tax_level', 'tax_prefix', 'formula', 'aldex_type','p_adj_method', 'prv_cut','abun_cut'),
                        Value = c(t, t2, formula, aldex_type,'fdr', prv_cut,abun_cut))
  
  return(list(results = aldex.res,input_ps = phyloseq,settings = settings))
}

#~~~~~~~~~

tidy_maaslin = function(phyloseq,formula,
                        abun_cut = 0.0001, prv_cut = 0.1, min_variance = 0.0,
                        normalization = "TSS", transform = "AST", analysis_method = "LM",
                        max_significance = 0.25, random_effects = NULL,
                        correction = "BH", standardize = F, cores = 1,
                        plot_heatmap = F, plot_scatter = F, heatmap_first_n = 50, reference = NULL){

  # phyloseq=ps;formula='TPN_control';
  # abun_cut = 0.001; prv_cut = 0.1; min_variance = 0.0;
  # normalization = "TSS"; transform = "AST"; analysis_method = "LM";
  # max_significance = 0.05; random_effects = NULL;
  # correction = "BH"; standardize = F; cores = 1;
  # plot_heatmap = F; plot_scatter = F; heatmap_first_n = 50; reference = NULL
  
  # Adapt formula for maas input
  fixed_effects = formula %>% str_split('[+]') %>% lapply(str_trim) %>% unlist()
  
  ## Taxonomic level
  t = which(!is.na(phyloseq@tax_table[1,])) %>% max # select lowest level
  t = colnames(phyloseq@tax_table)[t] # Name of tax level
  t2 = str_sub(t,1,1) %>% str_to_lower() # ex. 's' for Species
  
  o = phyloseq@otu_table %>% as.data.frame()
  s = phyloseq@sam_data %>% as.matrix() %>% as.data.frame()
  # Convert to numeric where applicable
  # Note: if the first value in a column is NA, then it won't convert properly.
  w = which(!is.na(as.numeric(s[1,]))) %>% suppressWarnings()
  s[, w] <- lapply(s[, w, drop = FALSE], function(x) as.numeric(x))
  
  # Run Maaslin2
  set.seed(421)
  maas.sp = Maaslin2(
    input_data = o, 
    input_metadata = s, 
    output = 'temp', 
    min_prevalence = prv_cut,
    max_significance = max_significance,
    normalization = normalization,
    transform = transform,
    analysis_method = analysis_method,
    fixed_effects = fixed_effects,
    correction = correction,
    standardize = standardize,
    plot_heatmap = plot_heatmap,
    plot_scatter = plot_scatter)
  
  # Removes Maaslin results folder
  unlink('temp',recursive = T) 
  
  # Format
  maas.res = format_maaslin_results(maas.sp,phyloseq,min_abundance = abun_cut, 
                                    min_prevalence = prv_cut,
                                    min_variance = min_variance, normalization = normalization,
                                    transform = transform)
  
  if(is.null(random_effects)){random_effects=NA}
  settings = data.frame(Parameter = c('tax_level', 'tax_prefix', 'formula', 'p_adj_method', 
                                      'prv_cut', 'abun_cut','min_variance',
                                      'normalization','transform','analysis_method',
                                      'max_significance','random_effects','standardize'),
                        Value = c(t, t2, formula, correction, 
                                  prv_cut,abun_cut,min_variance,
                                  normalization,transform,analysis_method,
                                  max_significance,random_effects,standardize))
  
  return(list(results=maas.res,input_ps = phyloseq,settings = settings))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# res is the results of running ancombc()
# ps is the phyloseq object you're using for the DA analysis
# Level letter = 'f' for family, etc
format_ancom_results = function(res,ps,min_abundance,tlevel){
  # res = ancom.sp;ps = phyloseq %>% tax_glom(tlevel);min_abundance = abun_cut;tlevel = tlevel
  
  # RESULTS FORMATTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # First we need to update the column names for each part of the results, because they're currently all the same.
  colnames(res$res$lfc) = paste(colnames(res$res$lfc),'__lfc',sep='')
  colnames(res$res$se) = paste(colnames(res$res$se),'__se',sep='')
  colnames(res$res$W) = paste(colnames(res$res$W),'__W',sep='')
  colnames(res$res$p_val) = paste(colnames(res$res$p_val),'__p_val',sep='')
  colnames(res$res$q_val) = paste(colnames(res$res$q_val),'__q_val',sep='')
  colnames(res$res$diff_abn) = paste(colnames(res$res$diff_abn),'__diff_abn',sep='')
  
  # First, we'll use lapply to apply a function to each item in a list.
  # For each item, the row names will be converted into a column called 'OTU'.
  # We then convert to a tibble (similar to a data frame, but with more flexibility)
  # Finally, we use reduce to collapse the list into a single table using full_join.
  # Full_join will combine the tables by any shared columns (OTU, in this case) without
  # accidentally removing any data.
  
  ## Abundance filter (prev filter built in to ancombc())
  ps2 <- transform_sample_counts(ps, function(x) x / sum(x))
  o = ps2@otu_table %>% as.matrix() %>% as.data.frame() %>% 
    dplyr::select(all_of(colnames(res$feature_table)))
  num_abun = apply(o, 1, function(x) mean(x,na.rm=T))
  
  tax_del = c(which(num_abun < min_abundance))
  if (length(tax_del) > 0){
    o.filt = o[-tax_del,]
  } else {
    o.filt = o
  }
  o.filt = ps2@tax_table %>% as.matrix %>% as.data.frame %>% filter(rownames(.) %in% rownames(o.filt))
  
  results = lapply(res$res,function(x) rownames_to_column(x,'OTU')) %>% 
    lapply(as.data.frame) %>% purrr::reduce(left_join) %>% as.matrix() %>% as_tibble() %>% 
    filter(OTU %in% o.filt[[tlevel]]) %>% # abundance filter
    dplyr::select(-contains('Intercept'))
  names(results)[which(names(results)=='OTU')] = tlevel
  
  # TAXONOMY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Add all taxonomic info - DEPRECATED WITH NEW ANCOMBC
  print('Fixing taxonomy')
  tax_info = update_taxonomy(ps)
  
  # Join taxonomy info to the results
  results.formatted = results %>% 
    left_join(tax_info[,c(1,ncol(tax_info))]) %>% # Get rid of other levels
    mutate(DA = 'ANCOM-BC') %>% 
    dplyr::select(-OTU) %>% 
    dplyr::select(names(tax_info)[ncol(tax_info)],DA,everything())
  
  # Convert to numeric where applicable
  w = which(!is.na(as.numeric(results.formatted[1,]))) %>% suppressWarnings()
  results.formatted[, w] <- lapply(results.formatted[, w, drop = FALSE], function(x) as.numeric(x))
  
  return(results.formatted)
}

#~~~~~~~~

# df.clr is the output of aldex.clr
# res.var isn't required for adlex.glm, but it'll add the string to the start of the stat column in order to match the other results tables. This'll help avoid issues when combining data if you used a type of aldex that doesn't add the variable name (ex. aldex.kw). ex. res.var = 'Fractionpos'
# min_prevalence is the minimum prevalence required to include the taxa. This should be
# the same value as was used to run aldex2.
format_aldex_results = function(res,ps,df.clr,res.var=NA,min_prevalence = 0.10,min_abundance=0.0001,
                                aldex_type = aldex_type){
  require(data.table)
  # res=aldex.sp;ps=ps;df.clr=x;
  # res.var=res.var; min_abundance=0.0001;
  # min_prevalence = prv_cut;aldex_type=aldex_type

  # First, recalculate the padj to only include taxa above the cutoff~~~~~~~~~~~
  ## Prevalence:
  s = ps@sam_data %>% as.matrix() %>% as.data.frame()
  o = ps@otu_table %>% as.matrix() %>% as.data.frame() %>% dplyr::select(all_of(names(df.clr@reads)))
  num_prop = apply(o, 1, function(x) sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))
  
  ## Abundance:
  ps2 <- transform_sample_counts(ps, function(x) x / sum(x))
  o = ps2@otu_table %>% as.matrix() %>% as.data.frame() %>% dplyr::select(all_of(names(df.clr@reads)))
  num_abun = apply(o, 1, function(x) mean(x,na.rm=T))
  
  tax_del = c(which(num_prop < min_prevalence), which(num_abun < min_abundance)) %>% unique() %>% sort
  if (length(tax_del) > 0){
    o.filt = o[-tax_del,]
  } else {
    o.filt = o
  }

  if('Row.names' %in% names(res)){
    res = res %>% column_to_rownames('Row.names')
  }
  
  res.filt = res[rownames(res) %in% rownames(o.filt),] %>% 
    dplyr::select(-contains('.BH'),-contains('Intercept'),-contains('.eBH'))
  srv = Vectorize(str_remove_all)
  names(res.filt) = srv(names(res.filt),'model.')
  #Test if output is kruskal-wallis stats
  if(aldex_type=='kw'){ # aldex.kw was used
    res.filt = res.filt %>% dplyr::select(-contains('glm'),-contains('effect[.]')) %>% 
      mutate(padj = p.adjust(res.filt$kw.ep,method='fdr')) %>% 
      dplyr::select(kw.ep,padj,everything())
    if(is.na(res.var)){
      names(res.filt)[1:2] = c('p_val','q_val')
    } else {
      names(res.filt)[1:2] = c('p_val','q_val')
      names(res.filt) = paste(res.var,names(res.filt),sep='__')
    }
  } else if(aldex_type=='spearman'){
    res.filt = res.filt %>% dplyr::select(contains('spearman')) %>% 
      rename(effect = spearman.erho,p_val = spearman.ep) %>% 
      mutate(q_val = p.adjust(p_val,method='fdr'))
    names(res.filt) = paste(res.var,names(res.filt),sep='__')
  } else if(aldex_type=='glm'){
    res.filt = res.filt %>% dplyr::select(-contains('win.'),-contains('overlap'))
    # First, separate by spaces and take the first value
    prefix = str_split(names(res.filt),' ')
    prefix = lapply(prefix,function(x) x = x[[1]][1]) %>% unlist()
    # Then separate by .
    prefix = str_split(prefix,'[.]')
    prefix = lapply(prefix,function(x) x = x[[1]][1]) %>% unlist()
    num_vars = length(unique(prefix))
    names(res.filt) = paste(prefix,c(
                            rep(c('estimate','se','t_stat','p_val'),num_vars),
                            rep(c('rab.all','diff.btw','diff.win','effect'),(ncol(res.filt)-(4*num_vars))/4)
                            ),sep='__')
    w = which(str_detect(names(res.filt),'p_val'))
    
    for(i in w){
      res.filt$padj = p.adjust(res.filt[,i],method = 'fdr')
      names(res.filt)[ncol(res.filt)] = paste(prefix[i],'q_val',sep='__')
    }
  } else {print('Error - wrong aldex type in format_aldex_results')}
  
  # TAXONOMY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Note that the taxa aren't listed as family names, but as a random OTU within that family. Let's add the family info and remove the OTUs.
  tax_info = update_taxonomy(ps)
  
  # Join taxonomy info to the results
  results.formatted = res.filt %>% rownames_to_column('OTU') %>% 
    left_join(tax_info[,c(1,ncol(tax_info))]) %>% # Get rid of other levels
    mutate(DA = 'ALDEx2') %>% 
    dplyr::select(-OTU) %>% 
    dplyr::select(names(tax_info)[ncol(tax_info)],DA,everything())
  
   return(results.formatted)
}

#~~~~~~~~~~

# res_loc is the location of the maaslin results - need it in order to load the transformed data.
# The other settings should be the same as was used for maaslin2 - we need them here in order to recalculate the input dataset that maaslin uses for statistical testing.
format_maaslin_results = function(res,ps,min_abundance = 0,
                                  min_prevalence = 0.10,min_variance = 0,
                                  normalization = 'TSS',transform='AST'){

  # res = maas.sp;ps=ps;min_abundance = 0.001;min_prevalence = 0.1;min_variance = 0;normalization = 'TSS';transform='AST'
  
  sig = res$results # Stats
  
  # TAXONOMY~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  tax_info = update_taxonomy(ps)
  
  # Maaslin changes any '-' to '.', so they need to be changed back if applicable.
  if(str_detect(sig$feature[1],'k__Bacteria[.]') & 
     str_detect(tax_info$OTU[1],'k__Bacteria-')){
    sig$feature = str_replace_all(sig$feature,'[.]','-')
  } else { # Otherwise, just fix names to match
    tax_info$OTU = make.names(tax_info$OTU)
  }
  
  # Fix features
  ## IF the features are OTUs, they often start with a number which will automatically trigger the prefix X to be added.
  ## Since we can't be completely sure that the OUT didn't start with an X by chance, we'll ADD X to the taxonomy table instead.
  # Anything where the first character is a number gets an X added
  tax_info$OTU[!is.na(as.numeric(str_sub(tax_info$OTU,1,1)))] = 
    paste('X',tax_info$OTU[!is.na(as.numeric(str_sub(tax_info$OTU,1,1)))],sep='') %>%
    suppressWarnings()
  
  tax_info = tax_info %>% dplyr::rename(feature = OTU)
  
  # Implement abundance filter
  psm = ps %>% microbiome::transform('compositional') %>% psmelt() %>% 
    group_by(OTU) %>% summarize(Abundance = mean(Abundance)) %>% 
    filter(Abundance>min_abundance) %>% 
    mutate(OTU = ifelse(!is.na(as.numeric(str_sub(OTU,1,1))),paste0('X',OTU),OTU))
  
  results.formatted = sig %>%
    filter(feature %in% psm$OTU) %>%# Abundance filter
    left_join(tax_info[,c(1,ncol(tax_info))]) %>% # Get rid of other levels
    dplyr::select(-feature) %>% 
    mutate(DA = 'MaAsLin2') %>% 
    dplyr::select(names(tax_info)[ncol(tax_info)],DA,everything()) %>% 
    dplyr::select(-c(metadata,value, N, N.not.zero)) %>% 
    dplyr::rename('beta' = coef, 'se'=stderr, 'p_val' = pval, 'q_val' = qval) %>% # To match other results
    pivot_longer(cols = c(beta,se,p_val,q_val),names_to = 'test',values_to = 'value') %>% 
    pivot_wider(names_from = c(name,test),values_from = value,names_sep = '__')
  
  return(results.formatted)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This function takes the results from the three DA tools and returns 5 things:
# 1. combined_res is the total results table. All equivalent columns are conserved. Beta/estimate/log fold change columns vary in meaning, so instead they are replaced with association_dir, which indicates the direction of the association.
# 2. significance_summary gives a table that indexes whether each bug is significant in a true/false manner. Total_sig gives the number of tools where the bug was significant, and sig_in_at_least_2 gives the final significant microbes, assuming that significant = q<alpha in at least 2 DA results (q = FDR-corrected p value).
# 3. sig_microbes outputs the bugs significant in at least 2 tools as a vector.
# 4. ps.clr is a clr-transformed phyloseq object which may be useful for plotting microbes.
# 5. settings gives the settings used for each tool.
# Inputs: the entire output of each formatting command is provided (not just the stats part). Alpha indicates the significance threshold for corrected p values, and stat specifies whether you want to consider qvals or pvals.
combine_da_results = function(ancom,aldex,maas,alpha=0.05,stat='q_val'){
  
  # ancom=ancom.sp;aldex=aldex.sp;maas=maas.sp;alpha=0.05;stat='q_val'
  
  results = full_join(ancom$results,aldex$results) %>%
    full_join(maas$results) %>% 
    dplyr::select(-contains(c('_W','diff_abn','t_stat'))) %>% 
    pivot_longer(cols = contains(c('lfc','se','p_val','q_val','beta','effect')),
                 names_to = c('Variable','Statistic'),
                 values_to = 'Value',
                 names_sep = '__') %>% 
    mutate(Statistic = replace(Statistic, 
                               Statistic %in% c('lfc','beta','effect'),'association_dir')) %>% 
    filter(!is.na(Value)) %>% 
    pivot_wider(names_from = Statistic, values_from = Value) %>% 
    mutate(association_dir = as.character(association_dir>0)) %>% 
    mutate(association_dir = replace(association_dir, association_dir=='TRUE','Positive'),
           association_dir = replace(association_dir, association_dir=='FALSE','Negative')) %>% 
    suppressWarnings() %>% suppressMessages()
  
  # Determine which microbes are significant, and in how many datasets
  if(stat=='q_val'){
    sig = results %>% mutate(sig = as.numeric(q_val<alpha))
  } else if(stat=='p_val'){
    sig = results %>% mutate(sig = as.numeric(p_val<alpha))
  } else {print('stat must be "q_val" or "p_val".')}
  sig = sig %>% 
    dplyr::select(all_of(names(results)[1]), DA,Variable,sig) %>% 
    pivot_wider(names_from = DA, values_from = sig)
  sig$total_sig = apply(sig,1,function(x)sum(as.numeric(x[c("ANCOM-BC","ALDEx2","MaAsLin2")]),na.rm=T))
  sig = sig %>% mutate(sig_in_at_least_2 = as.numeric(total_sig>=2))
  
  sig_microbes = sig[sig$sig_in_at_least_2==1,1] %>% unlist() %>% unname()
  sig_microbes = sig_microbes[!is.na(sig_microbes)]
  
  # Combine settings
  s1 = ancom$settings %>% rename('ANCOM-BC' = Value)
  s2 = aldex$settings %>% rename('ALDEx2' = Value)
  s3 = maas$settings %>% rename('MaAsLin2' = Value)
  S = full_join(s1,s2) %>% full_join(s3)
  
  # CLR-transformed phyloseq object
  require(microbiome)
  df = microbiome::transform(ancom$input_ps,'clr')

  return(list(combined_res = results, significance_summary = sig, subsetted_ps = sig_microbes, 
              ps.clr = df,settings=S))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This function takes the output from combine_res_results and produces a Venn diagram of the significant microbes across all 3 tools. 
# The filename is the name of the Venn diagram and should not contain a suffix (ex. '.png'). Feel free to add additional filepaths ex) if my working directory is C://Data and I want to save in C://Data/Analysis/Plots, I just need 'Analysis/Plots/filename'.
# If your combined_results summarizes the results of multiple variables, use var to specify which variable to plot.
da_venn_diagram = function(combined_results, var = NA, filename='da_venn_diagram', add_labels = T){
  require(VennDiagram)
  require(RColorBrewer)
  # combined_results=res; var = NA; filename='da_venn_diagram'
  
  temp = combined_results$significance_summary
  
  if(is.na(var)){
    if(length(unique(temp$Variable))>1){
      print('Warning: combined_results contains multiple variables, but no variable has been selected. The Venn diagram will combine the results from all variables.')
      print("Variables:")
      print(unique(temp$Variable))
    }
  } else {
    temp = temp %>% filter(Variable == var)
  }
  names(temp)[1] = 'test' # Since we don't know what level we're testing
  venn_data = list(`ANCOM-BC` = temp$test[!is.na(temp$`ANCOM-BC`) & temp$`ANCOM-BC`==1],
                   ALDEx2 = temp$test[!is.na(temp$ALDEx2) & temp$ALDEx2==1],
                   MaAsLin2 = temp$test[!is.na(temp$MaAsLin2) & temp$MaAsLin2==1])
  
  myCol <- brewer.pal(3, "Dark2")
  
  if(add_labels){
    venn.diagram(
    x = venn_data,
    category.names = c("ANCOM-BC" , "ALDEx2" , "MaAsLin2"),
    filename = paste(filename,'.png',sep=''),
    output=T,
    
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    alpha = 0.5,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 180),
    cat.dist = c(0.07, 0.07, 0.055),
    cat.fontfamily = "sans",
    rotation = 1
  )
  } else {
  venn.diagram(
    x = venn_data,
    category.names = c("ANCOM-BC" , "ALDEx2" , "MaAsLin2"),
    filename = paste(filename,'.png',sep=''),
    output=T,
    
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    alpha = 0.5,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0
  )
  }

  # The above command also produces an annoying log file each time it runs, so this automatically deletes it.
  unlink(paste(filename,'.png.',Sys.Date(),'*.log',sep=''))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This function takes a lot of inputs, but hey, it works!
# combined_results: the output of combine_da_results
# ps.var and res.var: Both of these describe your variable of interest. Since the DA tools add on one of the variable levels, you need to specify the original variable name (ps.var, ex. "Fraction") AND the variable that is output by the DA tools (res.var, ex. "Fractionpos").
# facet_var will facet wrap the plot (subdivide it based on a second variable)
# x_label: fix the x label if you want. Useful for continuous variables.
# save: do you want to save the plots?
# save_loc: if NA and save=T, plots will save in the working directory. otherwise, set as desired filepath relative to working directory.
plot_all_sig_microbes = function(combined_results,ps.var,res.var,facet_var='',xlabel = NA,stat='p_val',save=F,save_loc=NA,save_suffix=''){
  print('initialized')
  # combined_results=res2;ps.var='mindscore';res.var='mindscore';facet_var='';save=T;save_loc='Plots/Bugs';stat='p_val'

  # Find significant microbes
  temp = combined_results$significance_summary %>% 
    filter(Variable == res.var,
           sig_in_at_least_2==1)
  names(temp)[1] = 'temp'
  taxa_to_test = temp$temp
  
  level = combined_results$settings$`ANCOM-BC`[combined_results$settings$Parameter=='tax_prefix']
  
  # Add more taxonomy info where needed to resolve duplicate names (ex. two 'f__' present)
  tax_info = update_taxonomy(combined_results$ps.clr)
  combined_results$ps.clr@tax_table = tax_table(tax_info %>% column_to_rownames('OTU') %>% as.matrix())
  
  # Format the clr-transformed data
  temp = combined_results$ps.clr %>% psmelt() %>% 
    dplyr::select(Sample, all_of(names(combined_results$significance_summary)[1]),
                  all_of(ps.var),any_of(facet_var),Abundance)
  names(temp)[2] = 'taxa' # Second column is the taxonomic level of interest
  temp = temp %>% filter(taxa %in% taxa_to_test) %>% 
    pivot_wider(names_from=taxa,values_from=Abundance)
  
  # Decide whether plots will be linear or boxplots.
  if(is.character(temp[[ps.var]]) | is.factor(temp[[ps.var]])){
    type = 'boxplot'
  } else if(is.numeric(temp[[ps.var]])|is.integer(temp[[ps.var]])){
    type = 'linear'
  } else {
    print('Error: Structure of chosen variable is not supported. 
            Please change to integer, numeric, character, or factor.')
    break
  }
  
  # plotted_pvals
  
  plotted_pvals = combined_results$combined_res
  names(plotted_pvals)[1] = 'temp'
  plotted_pvals = plotted_pvals %>% 
    filter(temp %in% taxa_to_test,
           Variable == res.var)
  if(stat=='q_val'){
    plotted_pvals = plotted_pvals %>% dplyr::select(temp,DA,q_val)
    plotted_pvals$symb = sapply(plotted_pvals$q_val,
                                function(q){
                                  ifelse(q>0.05,'',
                                         ifelse(q>0.01,"*",
                                                ifelse(q>0.001,'**',
                                                       ifelse(q>0.0001,'***','****'))))
                              })
  } else if(stat=='p_val'){
    plotted_pvals = plotted_pvals %>% dplyr::select(temp,DA,p_val)
    plotted_pvals$symb = sapply(plotted_pvals$p_val,
                                function(p){
                                  ifelse(p>0.05,' ',
                                         ifelse(p>0.01,"*",
                                                ifelse(p>0.001,'**',
                                                       ifelse(p>0.0001,'***','****'))))
                                })
  } else {stop("'stat' should be q_val or p_val")}
  plotted_pvals$symb_pre = sapply(plotted_pvals$DA,
                                  function(d){
                                    ifelse(d=='ANCOM-BC','AN:',
                                           ifelse(d=='ALDEx2','AL:','MA:')) 
                                  })

  plotted_pvals = plotted_pvals %>%
    dplyr::select(temp,DA,symb,symb_pre) %>% 
    pivot_wider(names_from = DA,values_from = c(symb,symb_pre))
  plotted_pvals$label1 = paste(plotted_pvals$`symb_pre_ANCOM-BC`,
                              plotted_pvals$symb_pre_ALDEx2,
                              plotted_pvals$symb_pre_MaAsLin2,
                              sep='\n')
  plotted_pvals$label2 = paste(plotted_pvals$`symb_ANCOM-BC`,
                               plotted_pvals$symb_ALDEx2,
                               plotted_pvals$symb_MaAsLin2,
                               sep='\n')
  
  for(t in taxa_to_test){
    # t = taxa_to_test[1]
    if(type=='boxplot'){
      # t=taxa_to_test[1]
      
      names(temp)[names(temp)==ps.var] = 'Grouping'
      min = min(temp[[t]],na.rm=T)
      max = max(temp[[t]],na.rm=T)
      
      L1 = plotted_pvals$label1[plotted_pvals$temp==t]
      L2 = plotted_pvals$label2[plotted_pvals$temp==t]
      
      p = temp %>% 
        ggplot(aes(x=Grouping,y=.data[[t]],fill=Grouping)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(height = 0,width = 0.2) +
        theme_classic(base_size = 20) +
        xlab('') +
        scale_fill_discrete(name = ps.var) + 
        geom_text(x=1.45,y = max*1.3,label=L1,size = 6,lineheight=0.8,hjust=1) +
        geom_text(x=1.55,y = max*1.26,label=L2,size = 6,lineheight=0.8,hjust=0) +
        ylim(min*0.95,1.4*max) + 
        geom_segment(aes(x = 1, y = max*1.07, xend = 2, yend = max*1.07))+
        geom_segment(aes(x = 1, y = max*1.04, xend = 1, yend = max*1.07))+
        geom_segment(aes(x = 2, y = max*1.04, xend = 2, yend = max*1.07))
      
    } else {
      
      min = min(temp[[t]],na.rm=T)
      max = max(temp[[t]],na.rm=T)
      
      L1 = plotted_pvals$label1[plotted_pvals$temp==t]
      L2 = plotted_pvals$label2[plotted_pvals$temp==t]
      
      jitter_x = sd(temp[[ps.var]],na.rm=T)/10
      # Labels will be added around this x value:
      x_loc = min(temp[[ps.var]])+(max(temp[[ps.var]])-min(temp[[ps.var]]))/6
      
      p = temp %>% ggplot(aes(x=.data[[ps.var]],y=.data[[t]]),col='black') +
        geom_jitter(height = 0, width = jitter_x) +
        geom_smooth(method='lm',col='black') +
        ylim(min*0.95,1.35*max) + 
        theme_classic(base_size = 20) + 
        geom_text(x=x_loc-jitter_x,y = max*1.25,label=L1,size = 6,lineheight=0.8,hjust=1) +
        geom_text(x=x_loc+jitter_x,y = max*1.23,label=L2,size = 6,lineheight=0.8,hjust=0)
    }
    if(!is.na(xlabel)){
      p = p + xlab(xlabel)
    }
    if(facet_var != ''){
      p = p + facet_wrap(facet_var) + ylim(min*0.95,1.50*max)
      w = 8
      f = paste('_',facet_var,sep='')
    } else {
      w = 6
      f=''
    }
    if(save==T){
      if(!is.na(save_loc)){
        ggsave(paste(save_loc,'/',ps.var,'_',t,f,save_suffix,'.png',sep=''),plot = p, height = 5, width = w)
      } else {
        ggsave(paste(ps.var,'_',t,f,save_suffix,'.png',sep=''),plot = p, height = 5, width = w)
      }
      print(p)
    } else { 
      print(p)
    }
  }
}
