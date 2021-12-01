############################## Clean Up Environment, Load Libraries, and Load Custom Functions ##############################
# clean up environment
rm(list = ls())

# set working directory
# setwd("~/Dropbox/MBQC Manuscript/Code")

# load libraries
library(xgboost)
library(glmnet)
library(magrittr)
library(Matrix)
library(tidyr)
library(dplyr)
# library(RDS)
library(data.table)
library(batchtools)
library(ggplot2)
library(knitr)
library(gridExtra)
library(RColorBrewer)


# source custom functions
source("all_functions.R")

############################## Extract and Format Data Published in Nature Biotechnology ##############################
# load and format data posted on Nature Biotech site
if(!dir.exists("registry_for_mbqc_paper_data_extraction")){
  paper_extractor_reg <- makeRegistry(file.dir = "registry_for_mbqc_paper_data_extraction",
                                      packages = c("magrittr",
                                                   "glmnet",
                                                   "tidyr",
                                                   "dplyr",
                                                   # "RDS",
                                                   "readr",
                                                   "feather",
                                                   "R.utils",
                                                   "data.table"),
                                      source = "all_functions.R")

  batchMap(fun = paper_extractor,
           dummy_var = "hey there",
           reg = paper_extractor_reg)


  submitJobs(reg = paper_extractor_reg,
             resources = list(memory = "15G"))

  waitForJobs()

}

############################## Extract and Format Data Published at HMP MBQC site ##############################
# load and format data saved on hmp site
if(!dir.exists("registry_for_hmp_site_data_extraction")){
  site_extractor_reg <- makeRegistry(file.dir = "registry_for_hmp_site_data_extraction",
                                     packages = c("magrittr",
                                                  "glmnet",
                                                  "tidyr",
                                                  "dplyr",
                                                  "RDS",
                                                  "R.utils",
                                                  "data.table"),
                                     source = "all_functions.R")



  batchMap(fun = site_extractor,
           dummy_var = "hey there",
           reg = site_extractor_reg)


  submitJobs(reg = site_extractor_reg,
             resources = list(memory = "15G"))

  waitForJobs()

}


############################## Validate Data Extraction Workflow ##############################
# load metadata extracted from Nature Biotech via DSC's workflow (above)
# meta_paper_dsc <-
#   readRDS("mbqc_paper_data_products/meta_only") %>%
#   data.table
#
# # load metadata extracted from HMP MBQC site via DSC's workflow (also above)
# meta_site_dsc <-
#   readRDS("hmp_mbqc_site_data_products/meta_only") %>%
#   data.table
#
# # # load metadata extracted from HMP MBQC site via ADW's workflow (separate code)
# # load("meta_only.RData") #loads object 'meta_only'
# # meta_site_adw <-
# #   meta_only %>%
# #   data.table
# # rm(meta_only)
#
# ### Validation steps
#
# #1: check exact equality of 3 extracted data sources
# if(identical(mp_dsc_dim,ms_dsc_dim)){
#   print("Paper data (DSC extraction workflow) identical to HMP site data (DSC extraction workflow).")
# } else{
#   warning("Paper data (DSC extraction workflow) not identical to HMP site data (DSC extraction workflow).")
# }
#
# if(identical(mp_dsc_dim,ms_adw_dim)){
#   print("Paper data (DSC extraction workflow) identical to HMP site data (ADW extraction workflow).")
# } else{
#   warning("Paper data (DSC extraction workflow) not identical to HMP site data (ADW extraction workflow).")
# }
#
# if(identical(ms_dsc_dim,ms_adw_dim)){
#   print("HMP site data (DSC extraction workflow) identical to HMP site data (ADW extraction workflow).")
# } else{
#   warning("HMP site data (DSC extraction workflow) not identical to HMP site data (ADW extraction workflow).")
# }
#
# #2: check dimensions of data sources
# mp_dsc_dim <- meta_paper_dsc %>% dim
# ms_dsc_dim <- meta_site_dsc %>% dim
# ms_adw_dim <- meta_site_adw %>% dim
#
# #3: check variable names
# mp_dsc_vars <- meta_paper_dsc %>% colnames
# ms_dsc_vars <- meta_site_dsc %>% colnames
# ms_adw_vars <- meta_site_adw %>% colnames
#
# union_vars <- mp_dsc_vars %>% #union of all variable names
#   (function(x) union(x,ms_dsc_vars)) %>%
#   (function(x) union(x, ms_adw_vars))
#
# intersect_vars <- mp_dsc_vars %>% #intersection of all variable names
#   (function(x) intersect(x,ms_dsc_vars)) %>%
#   (function(x) intersect(x, ms_adw_vars))
#
# unshared_vars <- setdiff(union_vars,intersect_vars)
#
# if(length(unshared_vars) == 0){
#   print("All data sources contain same variables")
# } else{
#   warning(paste("The following variables are not shared between all data sources: \n",
#               paste(unshared_vars, collapse = " \n")))
# }
#
# ### Remove variables not shared between all sources
# mp_dsc_var_remove <- intersect(mp_dsc_vars, unshared_vars)
# meta_paper_dsc %<>%
#   select(-!!mp_dsc_var_remove)
#
# ms_dsc_var_remove <- intersect(ms_dsc_vars, unshared_vars)
# meta_site_dsc %<>%
#   select(-!!ms_dsc_var_remove)
#
# ms_adw_var_remove <- intersect(ms_adw_vars, unshared_vars)
# meta_site_adw %<>%
#   select(-!!ms_adw_var_remove)
#
# ### make column order identica
# meta_paper_dsc %<>% select(!!intersect_vars)
# meta_site_dsc %<>% select(!!intersect_vars)
# meta_site_adw %<>% select(!!intersect_vars)
#
# ### Remove bioinformatics lab 10, which is not mentioned in MBQC paper
# meta_paper_dsc %<>%
#   filter(dry_lab != "BL-10") %<>%
#   data.table
#
# ### Check equality of data tables after removal of BL-10 and non-matching variables
#
# same.unique.vals <- function(datalist, #determine whether variable in arbitary number of
#                              #data sources has same unique values across sources
#                          variable){
#   # get number of data sources in list
#   n_data <- length(datalist)
#   # create list of vectors of unique values
#   valuelist <- lapply(datalist, function(x)
#     x %>% select(!!variable) %>% unlist %>% unique)
#   # get all unique values (aggregating across data sources)
#   uni_values <- valuelist %>% unlist %>% unique
#   # get length of unique value vector
#   uni_length <- length(uni_values)
#   # get length of vector of unique values form each data source
#   length_vec <- sapply(valuelist, length)
#
#   if(max(length_vec) == min(length_vec)){ #if the max and the min of lengths are equal
#     if(max(length_vec)==uni_length){ #and equal to number of unique values
#       return(TRUE) # return true
#     } else{
#       return(FALSE) # otherwise return false
#     }
#   } else{
#     return(FALSE) # also return false if max != min
#   }
#
# }
#
# # determine which if any variables have discrepant unique values across data sources
# colnames(meta_paper_dsc) %>% #all column names are identical, so data source doesn't matter in this line
#   sapply(function(x) same.unique.vals(list(meta_paper_dsc,
#                       meta_site_dsc,
#                       meta_site_adw),
#                  x)) %>%
#   (function(x) !x) %>%
#   which()
#
# meta_site_adw %>%
#   with(table(blinded_lab, sequencing_wetlab, useNA = "ifany"))
#




############################## Apply Exclusion Criteria ##############################
# load metadata
metadata <- readRDS("mbqc_paper_data_products/meta_only")

wet_labs <- metadata$sequencing_wetlab %>% unique
wet_labs %<>% (function(x) x[!is.na(x)])

metadata %<>% filter(sequencing_wetlab %in% wet_labs) #filter out NA wet labs

# count number of extraction/sequencing protocols in non-pre-extracted samples
# from each sequencing wet lab
num_unique_protocols <- numeric(length(wet_labs))
for(i in 1:length(wet_labs)){
  num_unique_protocols[i] <- metadata %>%
    filter(pre.extracted_bool == "FALSE") %>%
    filter(sequencing_wetlab == wet_labs[i]) %>%
    filter(
      ((extraction_wetlab == wet_labs[i]))|
        (sequencing_wetlab == "HL-F" & extraction_wetlab == "HL-G")) %>%
    dplyr::select(extraction_kit_maker,
           extraction_kit_model,
           extraction_wetlab,
           homoginizer_maker,
           homoginizer_used,
           lab_transportation_method,
           seq_chem_version,
           seq_machine,
           storage_method_amplification_sequencing,
           storage_method_extraction_amplification,
           storage_method_receipt_processing,
           thaw_method_amplification,
           thaw_method_processing,
           thaw_method_sequencing,
           X16S_primer,
           X16S_regions
    ) %>%
    unique() %>%
    dim() %>%
    (function(x) x[1])
}

if(max(num_unique_protocols) == 1){
  print("Metadata indicates that no sequencing lab ran multiple protocols on non-pre-extracted samples")
} else{



  metadata$sequencing_wetlab_archived <- metadata$sequencing_wetlab

  for(w_l in wet_labs[num_unique_protocols>1]){
    length_unique <- metadata %>%
      filter(pre.extracted_bool == "FALSE") %>%
      filter(sequencing_wetlab == w_l) %>%
      filter((extraction_wetlab == w_l)|
               (sequencing_wetlab == "HL-F" & extraction_wetlab == "HL-G")) %>%
      dplyr::select(extraction_kit_maker,
             extraction_kit_model,
             extraction_wetlab,
             homoginizer_maker,
             homoginizer_used,
             lab_transportation_method,
             seq_chem_version,
             seq_machine,
             storage_method_amplification_sequencing,
             storage_method_extraction_amplification,
             storage_method_receipt_processing,
             thaw_method_amplification,
             thaw_method_processing,
             thaw_method_sequencing,
             X16S_primer,
             X16S_regions) %>%
      apply(2, function(x) length(unique(x)))

    key_pair <- metadata %>%
      filter(pre.extracted_bool == "FALSE") %>%
      filter(sequencing_wetlab == w_l) %>%
      filter((extraction_wetlab == w_l)|
               (sequencing_wetlab == "HL-F" & extraction_wetlab == "HL-G")) %>%
      dplyr::select(extraction_kit_maker,
             extraction_kit_model,
             extraction_wetlab,
             homoginizer_maker,
             homoginizer_used,
             lab_transportation_method,
             seq_chem_version,
             seq_machine,
             storage_method_amplification_sequencing,
             storage_method_extraction_amplification,
             storage_method_receipt_processing,
             thaw_method_amplification,
             thaw_method_processing,
             thaw_method_sequencing,
             X16S_primer,
             X16S_regions
      ) %>%
      dplyr::select(!!which(length_unique>1)[1]) %>% #this works bc max no. of protocols is 2
      ( function(x) list(colnames(x), unique(x)))


    suffix <-  metadata %>%
      filter(pre.extracted_bool == "FALSE") %>%
      filter(sequencing_wetlab == w_l) %>%
      filter((extraction_wetlab == w_l)|
               (sequencing_wetlab == "HL-F" & extraction_wetlab == "HL-G")) %>%
      dplyr::select(!!key_pair[[1]]) %>%
      apply(2, function(x) sapply(x, function(y) ifelse(y== key_pair[[2]][1,1],
                                                        1,2)))

    curr_wet_lab <- metadata %>%
      filter(pre.extracted_bool == "FALSE") %>%
      filter(sequencing_wetlab == w_l) %>%
      filter((extraction_wetlab == w_l)|
               (sequencing_wetlab == "HL-F" & extraction_wetlab == "HL-G")) %>%
      dplyr::select(sequencing_wetlab) %>%
      unlist()

    metadata[(metadata$pre.extracted_bool == "FALSE") &
               (metadata$sequencing_wetlab == w_l) &
               !is.na(metadata$sequencing_wetlab) &
               !is.na(metadata$pre.extracted_bool),"sequencing_wetlab"] <-
      paste(curr_wet_lab,suffix,sep="_")

    saveRDS(suffix, file = paste(w_l,"suffix",sep="_"))

  }
}
wet_labs <- metadata %>%
  filter(pre.extracted_bool == "FALSE") %>%
  dplyr::select(sequencing_wetlab) %>%
  unlist %>%
  unique

extraction_wetlabs <- metadata %>%
  filter(!is.na(extraction_wetlab)) %>%
  dplyr::select(extraction_wetlab) %>%
  unlist %>%
  unique

dry_labs <- metadata %>%
  filter(!is.na(dry_lab)) %>%
  dplyr::select(dry_lab) %>%
  unlist %>%
  unique

#calculate proportion of samples present (out of total run - should be multiple of 53 for each sequencing lab)
completeness_non_preextracted <- metadata %>%
  filter(sequencing_wetlab %in% wet_labs) %>% #include only samples from valid sequencing labs
  filter(extraction_wetlab %in% extraction_wetlabs) %>% #include only samples from valid extraction labs
  filter(pre.extracted_bool == FALSE) %>% #include only non-pre-extracted samples
  # filter(as.character(sequencing_wetlab_archived) == #include only samples extracted and sequenced
  #          as.character(extraction_wetlab)) %>% #at same lab
  with(table(dry_lab,sequencing_wetlab)) %>% #tabulate by dry lab and sequencing lab
  apply(c(1,2), function(x) if(x==0){return(0)} else{
    return((x/53)/ceiling(x/53))}) #get completeness of pre-extracted samples

#determine which combinations of sequencing and dry labs
#(for samples extracted at sequencing lab only)
#have all specimens present
metadata$sequencing_wetlab %<>% as.factor()
all_spec_non_preextracted <- metadata %>%
  filter(sequencing_wetlab %in% wet_labs) %>% #include only samples from valid sequencing labs
  filter(extraction_wetlab %in% extraction_wetlabs) %>% #include only samples from valid extraction labs
  filter(pre.extracted_bool == FALSE) %>% #include only non-pre-extracted samples
  # filter(as.character(sequencing_wetlab) == #include only samples extracted and sequenced
  #          as.character(extraction_wetlab)) %>% #at same lab
  group_by(sequencing_wetlab) %>%
  group_by(dry_lab,add = TRUE) %>%
  summarize(num_unique_spec = length(unique(specimen)))  %>% #get number of unique specimens by wet lab x dry lab
  ungroup %>%
  filter(num_unique_spec >= 22) %>%
  with(table(dry_lab,sequencing_wetlab))

non_redundant_cols <- which(!(colnames(all_spec_non_preextracted) %in% c("HL-F","HL-N")))

all_spec_non_preextracted <- all_spec_non_preextracted[,non_redundant_cols]

satisfactory <- all_spec_non_preextracted * apply(completeness_non_preextracted, c(1,2),
                                                  function(x) as.numeric(x>= .75))

### only consider wet labs that are satisfactory (have sufficiently complete data) for some dry lab
candidate_wet_labs <- satisfactory %>%
  apply(2,max) %>%
  (function(x) names(x[x>0]))

### update satisfactory
satisfactory <- satisfactory[,candidate_wet_labs]

### only consider dry labs that are satisfactory (have sufficiently complete data) for some wet lab
candidate_dry_labs <- satisfactory %>%
  apply(1, max) %>%
  (function(x) names(x[x>0]))

### function to determine whether a given combination of dry and wet labs
### meet exclusion criteria
is.satisfactory <- function(wet_labs, dry_labs,satisfaction_matrix){
  ### function to determine whether a given set of wet and dry labs are jointly satisfactory
  if(length(dry_labs)<4){ #must have at least 4 dry labs
    return(FALSE)
  }
  # if we have four dry labs then we want to make sure we have
  # complete enough data for every wet lab x dry lab combination
  satisfaction <- satisfaction_matrix[dry_labs,wet_labs] %>% min
  if(satisfaction ==0){
    return(FALSE)} else if(satisfaction == 1){
      return(TRUE)
    } else{
      return(NA)
    }
}


dry_lab_combns <- lapply(4:length(candidate_dry_labs), 
                         function(x) combn(candidate_dry_labs,x, simplify = FALSE)) %>%
  unlist(recursive = FALSE)
num_dl_combs <- length(dry_lab_combns)

unsatisfied <- TRUE
num_wet_labs <- length(candidate_wet_labs)
while(unsatisfied){
  print(paste("Exploring combinations of", num_wet_labs,"wet labs",sep = " "))
  wet_lab_combns <- combn(candidate_wet_labs,num_wet_labs, simplify  = FALSE)
  
  
  num_wl_combs <- length(wet_lab_combns)
  
  wl_results <- sapply(wet_lab_combns, function(x){
    sapply(dry_lab_combns, function(y){
      return(is.satisfactory(x,y, satisfactory))}
    )
  }
  )
  
  unsatisfied <- (sum(wl_results )== 0)
  
  print(paste("Found", sum(wl_results), "satisfactory wet lab x dry lab combination(s)",paste = " "))
  num_wet_labs <- num_wet_labs -1
  
}

### How many of these results include both HL-F_1 and HL-F_2 (two protocols, but only one wet lab)
# count number of times HL-F shows up in each combination
num_F <- wet_lab_combns %>% 
  sapply(function(x) as.numeric(("HL-F_1" %in% x) + ("HL-F_2" %in% x)))


### Get lengths of dry lab combinations
dry_lab_combn_lengths <- sapply(dry_lab_combns,length)
### filter results and wet lab combinations by combinations with a solution
dry_lab_unfound <- TRUE
dl_length <- length(candidate_dry_labs)

while(dry_lab_unfound){
  print(paste("Searching for solutions with", dl_length, "dry labs", sep = " "))
  #filtering by combinations with no more than one instance of HL-F
  num_sol <- wl_results[dry_lab_combn_lengths == dl_length,num_F<2] %>%
    sum()
  
  print(paste("Found", num_sol, "solutions with", dl_length, "dry labs",sep = " "))
  dry_lab_unfound <- (num_sol == 0)
  if(dry_lab_unfound){
    dl_length <- dl_length - 1}
  
}

#get indices of wet labs that meet criteria (with some group of dry labs)
wet_lab_indices <- wl_results[dry_lab_combn_lengths == dl_length,num_F<2] %>% apply(2, sum) %>% 
  (function(x) which(x ==1))

# get wet lab combinations to examine
suitable_wl_combns <- wet_lab_combns[num_F<2][wet_lab_indices]

# get bioinformatics lab combinations corresponding to above wet lab combinations
dry_lab_indices <- wl_results[,num_F<2][,wet_lab_indices] %>%
  apply(1,sum) %>%
  (function(x) which(x>0))

suitable_dl_combns <- dry_lab_combns[dry_lab_indices]

#pull out labs that differ between wet lab sets
diff_labs <- suitable_wl_combns %>%
  (function(x) setdiff(union(x[[1]],x[[2]]), intersect(x[[1]],x[[2]])))

#calculate mean completeness across dry labs to be used
mean_wl_1 <- completeness_non_preextracted[suitable_dl_combns[[1]],diff_labs[1]] %>% mean
mean_wl_2 <- completeness_non_preextracted[suitable_dl_combns[[2]],diff_labs[[2]]] %>% mean

# choose the combination containing the wet lab with higher mean completeness
combn_to_choose <- which.max(c(mean_wl_1,mean_wl_2))

included_wet_labs <- suitable_wl_combns[[combn_to_choose]]
included_dry_labs <- suitable_dl_combns[[combn_to_choose]]

# duplicate wet and dry labs for use in hl-c validation 
included_wet_labs_legacy <- included_wet_labs
included_dry_labs_legacy <- included_dry_labs


############################## Aggregate on Taxonomy and Split Out Test and Training Data ########################################
#set seed
train_seed <- 54332
if(!dir.exists("registry_for_data_aggregation")){
  agg_reg <- makeRegistry(file.dir = "registry_for_data_aggregation",
                          packages = c("magrittr",
                                       "glmnet",
                                       "tidyr",
                                       "dplyr",
                                       "RDS",
                                       "data.table"),
                          source = "all_functions.R")
  
  agg_params <- expand.grid(wet_lab = included_wet_labs,
                            level = c("phylum",
                                      "order",
                                      "class",
                                      "genus",
                                      "family",
                                      "species", 
                                      "OTU"))
  
  batchMap(fun = split_test,
           wet_lab = agg_params[,1],
           level = agg_params[,2],
           more.args = list(included_dry_labs = included_dry_labs,
                            seed = train_seed),
           reg = agg_reg)
  
  
  submitJobs(reg = agg_reg,
             resources = list(memory = "5G"))
  
  waitForJobs()
  
}

############################## Select and train classifiers ##############################################

### Parameters for elastic net (glmnet) classification
glmnet_params <- expand.grid(included_wet_labs,
                             alpha = seq(0,1, by = .05),
                             c("specimen","specimen_type_collapsed"),
                             c("OTU",
                               "phylum",
                               "order",
                               "class",
                               "genus",
                               "family",
                               "species"))
# Parameters for boosted tree (xgb) classifiers
xgb_params <- expand.grid(included_wet_labs,
                          eta = .1,
                          6,
                          subsample = c(0.5, 1.0),
                          colsample_bytree = seq(.25,1,by=.25),
                          c("specimen","specimen_type_collapsed"),
                          c("OTU",
                            "species",
                            "genus",
                            "family",
                            "order",
                            "class",
                            "phylum"))


### sensitivity analyses using data from dry labs individually
included_wet_labs_archive <- included_wet_labs
#exclude wet labs with only one aliquot set for sensitivity analysis
included_wet_labs <- included_wet_labs[!(included_wet_labs %in% 
                                         c("HL-I","HL-K","HL-L"))]
glmnet_params_archive <- glmnet_params
glmnet_params <- glmnet_params %>% dplyr::filter(!(Var1 %in% 
                                                     c("HL-I","HL-K","HL-L"))) %>%
  filter(Var3 == "specimen")

xgb_params_archive <- xgb_params
xgb_params <- xgb_params %>% dplyr::filter(!(Var1 %in% 
                                                     c("HL-I","HL-K","HL-L"))) %>%
  filter(Var6 == "specimen")
for(single_dry_lab in included_dry_labs){
  glmnet_cross_validate(included_wet_labs, transformation = "proportion_diff",
                        dry_lab_subset = single_dry_lab,
                        nfolds = 10)
  
  train_glmnet(glmnet_params,
               included_wet_labs,
               transformation = "proportion_diff",
               dry_lab_subset = single_dry_lab)
  
  xgb_cross_validate(included_wet_labs,
                     transformation = "proportion_diff",
                     dry_lab_subset = single_dry_lab,
                     resources = list(slots = 10,
                                      memory = "3G",
                                      queue = "students.q"))
  
  train_xgb(xgb_params, included_wet_labs,
            transformation = "proportion_diff",
            dry_lab_subset = single_dry_lab,
            resources = list(slots = 10,
                             memory = "3G",
                             queue = "students.q"))
  
  glmnet_cross_validate(included_wet_labs, transformation = "clr_pseudo_1_diff",
                        dry_lab_subset = single_dry_lab,
                        nfolds = 10)
  
  train_glmnet(glmnet_params,
               included_wet_labs,
               transformation = "clr_pseudo_1_diff",
               dry_lab_subset = single_dry_lab)
  
  xgb_cross_validate(included_wet_labs,
                     transformation = "clr_pseudo_1_diff",
                     dry_lab_subset = single_dry_lab,
                     resources = list(slots = 10,
                                      memory = "3G",
                                      queue = "students.q"))
  
  train_xgb(xgb_params, included_wet_labs,
            transformation = "clr_pseudo_1_diff",
            dry_lab_subset = single_dry_lab,
            resources = list(slots = 10,
                             memory = "3G",
                             queue = "students.q"))
}

included_wet_labs <- included_wet_labs_archive
glmnet_params <- glmnet_params_archive


glmnet_cross_validate(included_wet_labs, transformation = "clr_pseudo_1_diff_scale")

glmnet_cross_validate(included_wet_labs, transformation = "clr_pseudo_1_ecdf")

glmnet_cross_validate(included_wet_labs, transformation = "proportion_diff_scale")

glmnet_cross_validate(included_wet_labs, transformation = "proportion_ecdf")

train_glmnet(glmnet_params,
             included_wet_labs,
             transformation = "proportion_diff_scale")

train_glmnet(glmnet_params,
             included_wet_labs,
             transformation = "clr_pseudo_1_diff_scale")

train_glmnet(glmnet_params,
             included_wet_labs,
             transformation = "proportion_ecdf")

train_glmnet(glmnet_params,
             included_wet_labs,
             transformation = "clr_pseudo_1_ecdf")

glmnet_cross_validate(included_wet_labs, transformation = "proportion")

train_glmnet(glmnet_params,
included_wet_labs,
transformation = "proportion")


glmnet_cross_validate(included_wet_labs, transformation = "proportion_diff")

train_glmnet(glmnet_params,
included_wet_labs,
transformation = "proportion_diff")




glmnet_cross_validate(included_wet_labs, transformation = "presence")

train_glmnet(glmnet_params,
included_wet_labs,
transformation = "presence")


glmnet_cross_validate(included_wet_labs, transformation = "clr_pseudo_1")

train_glmnet(glmnet_params,
included_wet_labs,
transformation = "clr_pseudo_1")


glmnet_cross_validate(included_wet_labs, transformation = "clr_pseudo_1_diff")

train_glmnet(glmnet_params,
included_wet_labs,
transformation = "clr_pseudo_1_diff")




xgb_cross_validate(included_wet_labs,
                   transformation = "proportion_diff_scale",
                   resources = list(slots = 10,
                                    memory = "3G",
                                    queue = "students.q"))

train_xgb(xgb_params, included_wet_labs,
          transformation = "proportion_diff_scale",
          resources = list(slots = 10,
               memory = "3G",
               queue = "students.q"))


xgb_cross_validate(included_wet_labs,
                   transformation = "clr_pseudo_1_diff_scale",
                   resources = list(slots = 10,
                                    memory = "3G",
                                    queue = "students.q"))

train_xgb(xgb_params, included_wet_labs,
          transformation = "clr_pseudo_1_diff_scale",
          resources = list(slots = 10,
                           memory = "3G",
                           queue = "students.q"))


xgb_cross_validate(included_wet_labs,
                   transformation = "proportion",
                   resources = list(slots = 10,
                                    memory = "3G",
                                    queue = "students.q"))

train_xgb(xgb_params, included_wet_labs,
          transformation = "proportion",
          resources = list(slots = 10,
                           memory = "3G",
                           queue = "students.q"))



xgb_cross_validate(included_wet_labs,
                   transformation = "proportion_diff",
                   resources = list(slots = 10,
                                    memory = "3G",
                                    queue = "students.q"))

train_xgb(xgb_params, included_wet_labs,
          transformation = "proportion_diff",
          resources = list(slots = 10,
                           memory = "3G",
                           queue = "students.q"))



xgb_cross_validate(included_wet_labs,
                   transformation = "clr_pseudo_1",
                   resources = list(slots = 10,
                                    memory = "3G",
                                    queue = "students.q"))

train_xgb(xgb_params, included_wet_labs,
          transformation = "clr_pseudo_1",
          resources = list(slots = 10,
                           memory = "3G",
                           queue = "students.q"))



xgb_cross_validate(included_wet_labs,
                   transformation = "clr_pseudo_1_diff",
                   resources = list(slots = 10,
                                    memory = "3G",
                                    queue = "students.q"))

train_xgb(xgb_params, included_wet_labs,
          transformation = "clr_pseudo_1_diff",
          resources = list(slots = 10,
                           memory = "3G",
                           queue = "students.q"))





xgb_cross_validate(included_wet_labs,
                   transformation = "presence")

train_xgb(xgb_params, included_wet_labs,
          transformation = "presence")


############################## Collect Results ##############################


misclass_glmnet <- lapply(c("proportion_diff", #"proportion_diff_scale",
                            "clr_pseudo_1_diff", #"clr_pseudo_1_diff_scale",
                            "presence"),
                       function(x) get_misclass_results("glmnet",x))

misclass_xgb <- lapply(c("proportion_diff", #"proportion_diff_scale",
                         "clr_pseudo_1_diff", #"clr_pseudo_1_diff_scale",
                         "presence"),
                       function(x) get_misclass_results("xgb",x))


misclass_glmnet %<>% (function(x) do.call(rbind,x))

misclass_xgb %<>% (function(x) do.call(rbind,x))

misclass <- rbind(misclass_glmnet, misclass_xgb)

misclass$cross_lab <- (as.character(misclass$predictor_lab) ==
                         as.character(misclass$predictee_lab)) %>%
  sapply(function(x) ifelse(x,"Within Lab","Cross Lab"))
# #
misclass$class_var[misclass$class_var == "specimen"] <- "Specimen"
misclass$class_var[misclass$class_var == "specimen_type_collapsed"] <- "Specimen Type"
levels(misclass$class_var) <- c("Specimen", "Specimen Type")

misclass$predictee_lab %<>% as.factor()
misclass$predictor_lab %<>% as.factor()

levels(misclass$predictor_lab)[levels(misclass$predictor_lab)  == "HL-F_1"] <- "HL-F"
levels(misclass$predictee_lab)[levels(misclass$predictee_lab)  == "HL-F_1"] <- "HL-F"


############################## Collect Classification Results Quoted in Text ##############################

### Abstract numbers
misclass %>%
  filter(class_var == "Specimen") %>%
  group_by(cross_lab, transformation,taxon) %>%
  summarize(median_class = 1 - median(misclass)) %>%
  filter(transformation == "proportion_diff") %>%
  filter(taxon == "genus") %>%
  View()

#### Used in results

### 3.1

# specimen and specimen type median within lab
misclass %>%
  filter(transformation == "proportion_diff") %>%
  filter(cross_lab == "Within Lab") %>%
  group_by(class_var, classifier) %>%
  summarize(median_misclass = median(misclass),
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75))
# specimen and specimen type median within lab by taxon
misclass %>%
  filter(transformation == "proportion_diff") %>%
  filter(cross_lab == "Within Lab") %>%
  filter(class_var == "Specimen") %>%
  group_by(class_var, classifier, taxon) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75)) %>% View

### 3.1.1

# within and across lab for boosted tree and elastic net, averaged across taxa
misclass %>%
  filter(transformation == "proportion_diff") %>%
  filter(class_var == "Specimen") %>%
  group_by(classifier, cross_lab) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75))

# within and across by taxon (OTU and phylum)
misclass %>%
  filter(transformation == "proportion_diff") %>%
  filter(class_var == "Specimen") %>%
  filter(taxon %in% c("OTU","phylum")) %>%
  group_by(classifier, taxon, cross_lab) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75)) %>% View()

# within and across by taxon (species)
misclass %>%
  filter(transformation == "proportion_diff") %>%
  filter(class_var == "Specimen") %>%
  filter(taxon %in% c("species")) %>%
  group_by(classifier, taxon, cross_lab) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75)) %>% View()

### 3.1A

# specimen and specimen type median within lab
misclass %>%
  filter(transformation == "proportion_diff_scale") %>%
  filter(cross_lab == "Within Lab") %>%
  group_by(class_var, classifier) %>%
  summarize(median_misclass = median(misclass),
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75))
# specimen and specimen type median within lab by taxon
misclass %>%
  filter(transformation == "proportion_diff_scale") %>%
  filter(cross_lab == "Within Lab") %>%
  filter(class_var == "Specimen") %>%
  group_by(class_var, classifier, taxon) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75)) %>% View

### 3.1.1

# within and across lab for boosted tree and elastic net, averaged across taxa
misclass %>%
  filter(transformation == "proportion_diff_scale") %>%
  filter(class_var == "Specimen") %>%
  group_by(classifier, cross_lab) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75))

# within and across by taxon (OTU and phylum)
misclass %>%
  filter(transformation == "proportion_diff_scale") %>%
  filter(class_var == "Specimen") %>%
  filter(taxon %in% c("OTU","phylum")) %>%
  group_by(classifier, taxon, cross_lab) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75)) %>% View()

# within and across by taxon (species)
misclass %>%
  filter(transformation == "proportion_diff_scale") %>%
  filter(class_var == "Specimen") %>%
  filter(taxon %in% c("species")) %>%
  group_by(classifier, taxon, cross_lab) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75)) %>% View()

### 3.2

# specimen median within and across lab
misclass %>%
  filter(transformation == "clr_pseudo_1_diff") %>%
  filter(class_var == "Specimen") %>% 
  group_by(cross_lab, classifier) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75))

# within and across by taxon (OTU and phylum)
misclass %>%
  filter(transformation == "clr_pseudo_1_diff") %>%
  filter(class_var == "Specimen") %>%
  filter(taxon %in% c("OTU","phylum")) %>%
  group_by(classifier, taxon, cross_lab) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75)) %>% View()

#within and across by classifier at species level
misclass %>%
  filter(transformation == "clr_pseudo_1_diff") %>%
  filter(class_var == "Specimen") %>%
  filter(taxon %in% c("species")) %>%
  group_by(classifier, taxon, cross_lab) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75)) %>% View()

### 3.2A

# specimen median within and across lab
misclass %>%
  filter(transformation == "clr_pseudo_1_diff_scale") %>%
  filter(class_var == "Specimen") %>% 
  group_by(cross_lab, classifier) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75))

# within and across by taxon (OTU and phylum)
misclass %>%
  filter(transformation == "clr_pseudo_1_diff_scale") %>%
  filter(class_var == "Specimen") %>%
  filter(taxon %in% c("OTU","phylum")) %>%
  group_by(classifier, taxon, cross_lab) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75)) %>% View()

#within and across by classifier at species level
misclass %>%
  filter(transformation == "clr_pseudo_1_diff_scale") %>%
  filter(class_var == "Specimen") %>%
  filter(taxon %in% c("species")) %>%
  group_by(classifier, taxon, cross_lab) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75)) %>% View()



### 3.3

# within and across by taxon (OTU and phylum)
misclass %>%
  filter(transformation == "presence") %>%
  filter(class_var == "Specimen") %>%
  filter(taxon %in% c("OTU","phylum")) %>%
  group_by(classifier, taxon, cross_lab) %>%
  summarize(median_misclass = median(misclass), 
            q1_misclass = quantile(misclass,.25),
            q3_misclass = quantile(misclass,.75)) %>% View()

### Discussion 

### 4.2.2

misclass %>%
  filter(transformation == "clr_pseudo_1_diff") %>%
  filter(classifier == "Boosted Tree") %>%
  filter(class_var == "Specimen") %>%
  filter(predictee_lab %in% c("HL-F","HL-B")) %>%
  filter(predictor_lab %in% c("HL-F","HL-B")) %>%
  filter(taxon %in% c("species")) %>%
  View()


### conclusion numbers

# correct classification on species-level data, within and across laboratories
misclass %>%
  filter(taxon == "species") %>%
  filter(class_var == "Specimen") %>%
  group_by(cross_lab,transformation) %>%
  summarize(med_correct_class =  median(1 -misclass),
            q1_correct_class =  quantile(1 -misclass,.25),
            q3_correct_class =  quantile(1 -misclass,.75)
            )

misclass %>%
  filter(taxon == "phylum") %>%
  filter(class_var == "Specimen") %>%
  group_by(cross_lab,transformation) %>%
  filter(transformation == "clr_pseudo_1_diff") %>%
  summarize(med_correct_class =  median(1 -misclass),
            q1_correct_class =  quantile(1 -misclass,.25),
            q3_correct_class =  quantile(1 -misclass,.75)
  )

misclass %>%
  filter(class_var == "Specimen") %>%
  filter(transformation == "clr_pseudo_1_diff") %>%
  group_by(cross_lab,transformation) %>%
  filter(taxon == "species") %>%
  summarize(med_correct_class =  median(1 -misclass),
            q1_correct_class =  quantile(1 -misclass,.25),
            q3_correct_class =  quantile(1 -misclass,.75)
  )


############################## Create Results and Discussion Figures ##############################

####### Results Figures

##### Proportion Difference Figures

misclass$predictor_lab %<>%
  as.character() %<>%
  factor(levels =
           as.character(misclass %>%
                          # filter(taxon == "phylum") %>%
                          filter(class_var == "Specimen") %>%
                          filter(transformation == "proportion_diff") %>%
                          filter(cross_lab == "Within Lab") %>%
                          group_by(predictor_lab) %>%
                          summarize(mean_within = mean(misclass)) %>%
                          with(predictor_lab[order(mean_within)]) %>%
                          as.character()
           )
  )

pdf(width = 4, height = 2, file = "proportion_within.pdf")
misclass %>%
  filter(transformation == "proportion_diff") %>%
  filter(cross_lab == "Within Lab") %>%
  # filter(classifier == "Boosted Tree") %>%
  ggplot(aes(x = taxon, y = misclass,
             color = predictee_lab,
             group = predictee_lab)) +
  geom_point(size = .5) +
  geom_line(size = .3) +
  facet_grid(classifier~class_var, scales = "free_y") +
  theme_bw(base_size = 8) +
  scale_color_brewer(palette = "Dark2") +
  guides(color = guide_legend(title = "Handling Laboratory",nrow = 4, byrow = 4, title.position="top")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        # legend.position = "bottom",
        plot.title = element_text(size = 8)
      
  ) +
  ylab("Within-Laboratory \nMisclassification") +
  xlab("Taxon")

dev.off()

misclass$transformation %<>% factor(levels = c("proportion_diff",
                                                  "clr_pseudo_1_diff",
                                                  "presence"))
misclass$classifier %<>% factor(levels = c("Elastic Net",
                                          "Boosted Tree"))

results.names <- c('Boosted Tree.proportion_diff' = "Proportion (Centered) \nBoosted Tree",
                   'Elastic Net.proportion_diff' = "Proportion (Centered) \nElastic Net",
                   'Boosted Tree.clr_pseudo_1_diff' = "CLR (Centered) \nBoosted Tree",
                   'Elastic Net.clr_pseudo_1_diff' = "CLR (Centered) \nElastic Net",
                   'OTU' = "OTU",
                   'genus' = "Genus",
                   'order' = "Order",
                   'phylum' = "Phylum")
                   


pdf(width = 6, height = 4.5, file = "results_all.pdf")
ggplot() +
  geom_point(data = misclass %>%
               filter(transformation %in% c("proportion_diff","clr_pseudo_1_diff")) %>%
               filter(class_var == "Specimen") %>%
               filter(taxon %in% c("OTU","genus","order","phylum")) %>%
               # filter(classifier == "Boosted Tree") %>%
               # filter(taxon == "genus") %>%
               # filter(cross_lab == "Within Lab") %>%
               mutate(grouper = paste(classifier,predictee_lab,
                                      sep = "_",
                                      collapse = "_")),
             aes(x= predictor_lab,
                 y = misclass,
                 color = predictee_lab,
                 shape = cross_lab,
                 size = cross_lab,
                 alpha = cross_lab),
             size = .5) +
  geom_line(data = misclass %>%
              filter(transformation %in% c("proportion_diff","clr_pseudo_1_diff")) %>%
              filter(class_var == "Specimen") %>%
              # filter(taxon == "OTU") %>%
              # filter(classifier  == "Boosted Tree") %>%
              filter(taxon %in% c("OTU","genus","order","phylum")) %>%
              filter(predictor_lab != predictee_lab) %>%
              mutate(grouper = paste(predictee_lab,
                                     sep = "_",
                                     collapse = "_")),
            aes(x = predictor_lab,
                y = misclass,
                color = predictee_lab,
                group = predictee_lab,
                alpha = cross_lab),
            size= .3
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_size_discrete(range = c(1,2)) +
  scale_alpha_discrete(range = c(.75,1)) +
  coord_cartesian(ylim= c(0,1)) +
  ylab("Misclassification") +
  xlab("Training Laboratory") +
  facet_grid(interaction(classifier,transformation) ~ taxon,
             labeller = as_labeller(results.names)) +
  # ggtitle("Proportion-Scale Specimen Classifier Performance by Taxon") +
  guides(color = guide_legend(title = "Test Laboratory", nrow = 4, byrow = 4, title.position = "top"),
         shape = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
                              override.aes = list(linetype  = 0)),
         size = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
                             override.aes = list(linetype  = 0)),
         alpha = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
                              override.aes = list(linetype  = 0))) +
  theme_bw(base_size = 8) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 8))
dev.off()

##### Proportion Difference-Scale Figures

misclass$predictor_lab %<>%
  as.character() %<>%
  factor(levels =
           as.character(misclass %>%
                          # filter(taxon == "phylum") %>%
                          filter(class_var == "Specimen") %>%
                          filter(transformation == "proportion_diff_scale") %>%
                          filter(cross_lab == "Within Lab") %>%
                          group_by(predictor_lab) %>%
                          summarize(mean_within = mean(misclass)) %>%
                          with(predictor_lab[order(mean_within)]) %>%
                          as.character()
           )
  )

pdf(width = 4, height = 2, file = "proportion_scale_within.pdf")
misclass %>%
  filter(transformation == "proportion_diff_scale") %>%
  filter(cross_lab == "Within Lab") %>%
  # filter(classifier == "Boosted Tree") %>%
  ggplot(aes(x = taxon, y = misclass,
             color = predictee_lab,
             group = predictee_lab)) +
  geom_point(size = .5) +
  geom_line(size = .3) +
  facet_grid(classifier~class_var, scales = "free_y") +
  theme_bw(base_size = 8) +
  scale_color_brewer(palette = "Dark2") +
  guides(color = guide_legend(title = "Handling Laboratory",nrow = 4, byrow = 4, title.position="top")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        # legend.position = "bottom",
        plot.title = element_text(size = 8)
        
  ) +
  ylab("Within-Laboratory \nMisclassification") +
  xlab("Taxon")

dev.off()

# pdf(width = 6, height = 3, file = "proportion_scale_all.pdf")
# ggplot() +
#   geom_point(data = misclass %>%
#                filter(transformation == c("proportion_diff_scale")) %>%
#                filter(class_var == "Specimen") %>%
#                filter(taxon %in% c("OTU","genus","order","phylum")) %>%
#                # filter(classifier == "Boosted Tree") %>%
#                # filter(taxon == "genus") %>%
#                # filter(cross_lab == "Within Lab") %>%
#                mutate(grouper = paste(classifier,predictee_lab,
#                                       sep = "_",
#                                       collapse = "_")),
#              aes(x= predictor_lab,
#                  y = misclass,
#                  color = predictee_lab,
#                  shape = cross_lab,
#                  size = cross_lab,
#                  alpha = cross_lab),
#              size = .5) +
#   geom_line(data = misclass %>%
#               filter(transformation == c("proportion_diff_scale")) %>%
#               filter(class_var == "Specimen") %>%
#               # filter(taxon == "OTU") %>%
#               # filter(classifier  == "Boosted Tree") %>%
#               filter(taxon %in% c("OTU","genus","order","phylum")) %>%
#               filter(as.character(predictor_lab) != as.character(predictee_lab)) %>%
#               mutate(grouper = paste(predictee_lab,
#                                      sep = "_",
#                                      collapse = "_")),
#             aes(x = predictor_lab,
#                 y = misclass,
#                 color = predictee_lab,
#                 group = predictee_lab,
#                 alpha = cross_lab),
#             size= .3
#   ) +
#   scale_color_brewer(palette = "Dark2") +
#   scale_size_discrete(range = c(1,2)) +
#   scale_alpha_discrete(range = c(.75,1)) +
#   coord_cartesian(ylim= c(0,1)) +
#   ylab("Misclassification") +
#   xlab("Predictor Laboratory") +
#   facet_grid(classifier ~ taxon) +
#   # ggtitle("Proportion-Scale Specimen Classifier Performance by Taxon") +
#   guides(color = guide_legend(title = "Predictee Laboratory", nrow = 4, byrow = 4, title.position = "top"),
#          shape = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
#                               override.aes = list(linetype  = 0)),
#          size = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
#                              override.aes = list(linetype  = 0)),
#          alpha = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
#                               override.aes = list(linetype  = 0))) +
#   theme_bw(base_size = 8) +
#   theme(axis.text.x=element_text(angle = 45, hjust = 1),
#         plot.title = element_text(size = 8))
# dev.off()

##### CLR Difference Figures

misclass$predictor_lab %<>%
  as.character() %<>%
  factor(levels =
           as.character(misclass %>%
                          # filter(taxon == "phylum") %>%
                          filter(class_var == "Specimen") %>%
                          filter(transformation == "clr_pseudo_1_diff") %>%
                          filter(cross_lab == "Within Lab") %>%
                          group_by(predictor_lab) %>%
                          summarize(mean_within = mean(misclass)) %>%
                          with(predictor_lab[order(mean_within)]) %>%
                          as.character()
           )
  )

color_tab <- color_tab[order(levels(misclass$predictor_lab), decreasing = FALSE),]

# 
# pdf(width = 6, height = 3, file = "clr_all.pdf")
# ggplot() +
#   geom_point(data = misclass %>%
#                filter(transformation == c("clr_pseudo_1_diff")) %>%
#                filter(class_var == "Specimen") %>%
#                filter(taxon %in% c("OTU","genus","order","phylum")) %>%
#                # filter(classifier == "Boosted Tree") %>%
#                # filter(taxon == "genus") %>%
#                # filter(cross_lab == "Within Lab") %>%
#                mutate(grouper = paste(classifier,predictee_lab,
#                                       sep = "_",
#                                       collapse = "_")),
#              aes(x= predictor_lab,
#                  y = misclass,
#                  color = predictee_lab,
#                  shape = cross_lab,
#                  size = cross_lab,
#                  alpha = cross_lab),
#              size = .5) +
#   geom_line(data = misclass %>%
#               filter(transformation == c("clr_pseudo_1_diff")) %>%
#               filter(class_var == "Specimen") %>%
#               # filter(taxon == "OTU") %>%
#               # filter(classifier  == "Boosted Tree") %>%
#               filter(taxon %in% c("OTU","genus","order","phylum")) %>%
#               filter(predictor_lab != predictee_lab) %>%
#               mutate(grouper = paste(predictee_lab,
#                                      sep = "_",
#                                      collapse = "_")),
#             aes(x = predictor_lab,
#                 y = misclass,
#                 color = predictee_lab,
#                 group = predictee_lab,
#                 alpha = cross_lab),
#             size = .3
#   ) +
#   scale_color_brewer(palette = "Dark2") +
#   scale_size_discrete(range = c(1,2)) +
#   scale_alpha_discrete(range = c(.75,1)) +
#   coord_cartesian(ylim= c(0,1)) +
#   ylab("Misclassification") +
#   xlab("Predictor Laboratory") +  
#   # ggtitle("Log-Ratio-Scale Specimen Classifier Performance by Taxon") +
#   facet_grid(classifier ~ taxon) +
#   guides(color = guide_legend(title = "Predictee Laboratory", nrow = 4, byrow = 4, title.position = "top"),
#          shape = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
#                               override.aes = list(linetype  = 0)),
#          size = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
#                              override.aes = list(linetype  = 0)),
#          alpha = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
#                               override.aes = list(linetype  = 0)))  +
#   theme_bw(base_size = 8) +
#   theme(axis.text.x=element_text(angle = 45, hjust = 1),
#         plot.title = element_text(size = 8))
# dev.off()

##### CLR Difference-Scale Figures

misclass$predictor_lab %<>%
  as.character() %<>%
  factor(levels =
           as.character(misclass %>%
                          # filter(taxon == "phylum") %>%
                          filter(class_var == "Specimen") %>%
                          filter(transformation == "clr_pseudo_1_diff_scale") %>%
                          filter(cross_lab == "Within Lab") %>%
                          group_by(predictor_lab) %>%
                          summarize(mean_within = mean(misclass)) %>%
                          with(predictor_lab[order(mean_within)]) %>%
                          as.character()
           )
  )

color_tab <- color_tab[order(levels(misclass$predictor_lab), decreasing = FALSE),]

# 
# pdf(width = 6, height = 3, file = "clr_scale_all.pdf")
# ggplot() +
#   geom_point(data = misclass %>%
#                filter(transformation == c("clr_pseudo_1_diff_scale")) %>%
#                filter(class_var == "Specimen") %>%
#                filter(taxon %in% c("OTU","genus","order","phylum")) %>%
#                # filter(classifier == "Boosted Tree") %>%
#                # filter(taxon == "genus") %>%
#                # filter(cross_lab == "Within Lab") %>%
#                mutate(grouper = paste(classifier,predictee_lab,
#                                       sep = "_",
#                                       collapse = "_")),
#              aes(x= predictor_lab,
#                  y = misclass,
#                  color = predictee_lab,
#                  shape = cross_lab,
#                  size = cross_lab,
#                  alpha = cross_lab),
#              size = .5) +
#   geom_line(data = misclass %>%
#               filter(transformation == c("clr_pseudo_1_diff_scale")) %>%
#               filter(class_var == "Specimen") %>%
#               # filter(taxon == "OTU") %>%
#               # filter(classifier  == "Boosted Tree") %>%
#               filter(taxon %in% c("OTU","genus","order","phylum")) %>%
#               filter(as.character(predictor_lab) != as.character(predictee_lab)) %>%
#               mutate(grouper = paste(predictee_lab,
#                                      sep = "_",
#                                      collapse = "_")),
#             aes(x = predictor_lab,
#                 y = misclass,
#                 color = predictee_lab,
#                 group = predictee_lab,
#                 alpha = cross_lab),
#             size = .3
#   ) +
#   scale_color_brewer(palette = "Dark2") +
#   scale_size_discrete(range = c(1,2)) +
#   scale_alpha_discrete(range = c(.75,1)) +
#   coord_cartesian(ylim= c(0,1)) +
#   ylab("Misclassification") +
#   xlab("Predictor Laboratory") +  
#   # ggtitle("Log-Ratio-Scale Specimen Classifier Performance by Taxon") +
#   facet_grid(classifier ~ taxon) +
#   # facet_grid(. ~ taxon) +
#   guides(color = guide_legend(title = "Predictee Laboratory", nrow = 4, byrow = 4, title.position = "top"),
#          shape = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
#                               override.aes = list(linetype  = 0)),
#          size = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
#                              override.aes = list(linetype  = 0)),
#          alpha = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
#                               override.aes = list(linetype  = 0)))  +
#   theme_bw(base_size = 8) +
#   theme(axis.text.x=element_text(angle = 45, hjust = 1),
#         plot.title = element_text(size = 8))
# dev.off()

##### Presence-Absence Figures

misclass$predictor_lab %<>%
  as.character() %<>%
  factor(levels =
           as.character(misclass %>%
                          filter(taxon == "phylum") %>%
                          filter(class_var == "Specimen") %>%
                          filter(transformation == "presence") %>%
                          # filter(classifier == "Boosted Tree") %>%
                          filter(cross_lab == "Within Lab") %>%
                          group_by(predictor_lab) %>%
                          summarize(mean_within = mean(misclass)) %>%
                          with(predictor_lab[order(mean_within)]) %>%
                          as.character()
           )
  )

pdf(width = 6, height = 3, file = "presence_all.pdf")
ggplot() +
  geom_point(data = misclass %>%
               filter(transformation == c("presence")) %>%
               filter(class_var == "Specimen") %>%
               filter(taxon %in% c("OTU","genus","order","phylum")) %>%
               # filter(classifier == "Boosted Tree") %>%
               # filter(taxon == "genus") %>%
               # filter(cross_lab == "Within Lab") %>%
               mutate(grouper = paste(classifier,predictee_lab,
                                      sep = "_",
                                      collapse = "_")),
             aes(x= predictor_lab,
                 y = misclass,
                 color = predictee_lab,
                 shape = cross_lab,
                 size = cross_lab,
                 alpha = cross_lab),
             size = .5) +
  geom_line(data = misclass %>%
              filter(transformation == c("presence")) %>%
              filter(class_var == "Specimen") %>%
              # filter(taxon == "OTU") %>%
              # filter(classifier  == "Boosted Tree") %>%
              filter(taxon %in% c("OTU","genus","order","phylum")) %>%
              filter(predictor_lab != predictee_lab) %>%
              mutate(grouper = paste(predictee_lab,
                                     sep = "_",
                                     collapse = "_")),
            aes(x = predictor_lab,
                y = misclass,
                color = predictee_lab,
                group = predictee_lab,
                alpha = cross_lab),
            size = .3
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_size_discrete(range = c(1,2)) +
  scale_alpha_discrete(range = c(.75,1)) +
  coord_cartesian(ylim= c(0,1)) +
  ylab("Misclassification") +
  xlab("Training Laboratory") +  
  facet_grid(classifier ~ taxon) +
  # ggtitle("Presence-Scale Specimen Classifier Performance by Taxon") +
  guides(color = guide_legend(title = "Test Laboratory", nrow = 4, byrow = 4, title.position = "top"),
         shape = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
                              override.aes = list(linetype  = 0)),
         size = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
                             override.aes = list(linetype  = 0)),
         alpha = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
                              override.aes = list(linetype  = 0))) +
  theme_bw(base_size = 8) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 8))
dev.off()


####### Discussion Figures

misclass_pretty <- misclass

misclass_pretty$transformation %<>% sapply(function(x) ifelse(x == "clr_pseudo_1_diff","CLR (Centered)",
                                                              ifelse(x == "clr_pseudo_1_diff_scale","CLR (Centered, Scaled)",
                                                              ifelse(x == "proportion_diff","Proportion (Centered)",
                                                              ifelse(x == "proportion_diff_scale", "Proportion (Centered, Scaled)",
                                                              ifelse(x == "proportion_ecdf","Proportion (Empirical CDF)",
                                                              ifelse(x == "clr_pseudo_1_ecdf","CLR (Empirical CDF)",
                                                                     "Presence-Absence")))))))

misclass_pretty$transformation  <- factor(misclass_pretty$transformation,
                                          levels = c("Proportion (Centered)",
                                                     "Proportion (Centered, Scaled)",
                                                     "Proportion (Empirical CDF)",
                                                     "CLR (Centered)",
                                                     "CLR (Centered, Scaled)",
                                                     "CLR (Empirical CDF)",
                                                     "Presence-Absence"))



misclass_pretty$grouper <- 
  apply(misclass_pretty,1, function(x) paste(x["cross_lab"],x["classifier"], sep = "_", collapse = "_"))

pdf(width = 6, height = 3, file = "median_misclass.pdf")
misclass_pretty %>%
  filter(transformation %in% c("Proportion (Centered)","CLR (Centered)")) %>%
  group_by %>%
  group_by(transformation, class_var, taxon, cross_lab, classifier) %>%
  summarize(median_misclass = median(misclass),
               grouper = unique(grouper))%>%

  ggplot() + 
  # geom_ribbon(aes(x = taxon, ymin = lower, ymax = upper, group = grouper, color = cross_lab,
                  # fill = cross_lab), alpha = .5)+
  # geom_jitter(aes(x = taxon, y = misclass, color = cross_lab, shape = classifier), alpha = .2,
  #             height = 0, width = .3,
  # data = misclass_pretty) +
  geom_line(aes(x = taxon, y = misclass, color = cross_lab, linetype = classifier,
                group = interaction(predictee_lab, predictor_lab, classifier)), alpha = .2, size = .2,
              data = misclass_pretty %>%
              filter(transformation %in% c("Proportion (Centered)","CLR (Centered)"))) +
  geom_line(aes(x = taxon, y = median_misclass, group = interaction(cross_lab, classifier), color = cross_lab,
                linetype = classifier),
            alpha = 1,
            size = .5) + 
  facet_grid(class_var ~ transformation) + 
  theme_bw(base_size = 8) + 
  # ggtitle("Classification Results by Transformation, Taxon, and Classification Task") +
  ylab("Misclassification Error") + 
  xlab("Taxon") + 
  scale_shape_manual(values = c(4,1))+
  theme(axis.text.x=element_text(angle = 45, hjust = 1)
        ) +
  guides(color = guide_legend(title = "Type of Prediction",nrow = 2, byrow = 2, title.position="top"),
         linetype = guide_legend(title = "Classifier",nrow = 2, byrow = 2, title.position="top"),
         fill = guide_legend(title = "Type of Prediction",nrow = 2, byrow = 2, title.position="top"),
         shape = guide_legend(title = "Classifier",nrow = 2, byrow = 2, title.position="top"))
dev.off()

levels(misclass_pretty$transformation) <- c("Proportion", "CLR", "Presence")

pdf(width = 6, height = 3, file = "best_of_misclass.pdf")
misclass_pretty %>%
  filter(predictor_lab %in% c("HL-J","HL-B","HL-H","HL-F")) %>%
  filter(predictee_lab %in% c("HL-J","HL-B","HL-H", "HL-F")) %>%
  filter(class_var == "Specimen") %>%
  filter(classifier == "Boosted Tree") %>%
  ggplot() +
  # ggtitle("Classification Results in Laboratories with Low Within-Laboratory Misclassification") +
  geom_line(aes(x = taxon, y = misclass, group = interaction(predictor_lab, predictee_lab, classifier), color = predictee_lab,
            linetype = cross_lab), size = .3) + 
  facet_grid(transformation~predictor_lab) + 
  theme_bw(base_size = 8) + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1))+
  scale_color_brewer(palette = "Dark2") +
  ylab("Misclassification Error") +
  xlab("Taxon") +
  guides(linetype = guide_legend(title = "Type of Prediction",nrow = 2, byrow = 2, title.position="top"),
         color = guide_legend(title = "Test Laboratory",nrow = 2, byrow = 2, title.position="top")) 
dev.off()


############################## Supplementary Figures ##############################

### Load results including uncentered proportion and CLR
misclass_glmnet <- lapply(c("proportion",
                            "proportion_diff",
                            "clr_pseudo_1",
                            "clr_pseudo_1_diff",
                            "presence"),
                          function(x) get_misclass_results("glmnet",x))

misclass_xgb <- lapply(c("proportion",
                         "proportion_diff",
                         "clr_pseudo_1",
                         "clr_pseudo_1_diff",
                         "presence"),
                       function(x) get_misclass_results("xgb",x))


misclass_glmnet %<>% (function(x) do.call(rbind,x))

misclass_xgb %<>% (function(x) do.call(rbind,x))

misclass <- rbind(misclass_glmnet, misclass_xgb)

misclass$cross_lab <- (as.character(misclass$predictor_lab) ==
                         as.character(misclass$predictee_lab)) %>%
  sapply(function(x) ifelse(x,"Within Lab","Cross Lab"))
# #
levels(misclass$class_var) <- c("Specimen", "Specimen Type")

levels(misclass$predictor_lab)[levels(misclass$predictor_lab)  == "HL-F_1"] <- "HL-F"
levels(misclass$predictee_lab)[levels(misclass$predictee_lab)  == "HL-F_1"] <- "HL-F"


make_supp_results_fig <- function(misclass_data,
                          trans,
                          class_type){
  
misclass_data$predictor_lab %<>%
    as.character() %<>%
    factor(levels =
             as.character(misclass %>%
                            filter(class_var == class_type) %>%
                            filter(transformation == trans) %>%
                            filter(cross_lab == "Within Lab") %>%
                            group_by(predictor_lab) %>%
                            summarize(mean_within = mean(misclass)) %>%
                            with(predictor_lab[order(mean_within)]) %>%
                            as.character()
             )
    )
  

  gg <- ggplot() +
    geom_point(data = misclass_data %>%
                 filter(transformation == trans) %>%
                 filter(class_var == class_type) %>%
                 # filter(taxon %in% c("OTU","genus","order","phylum")) %>%
                 # filter(classifier == "Boosted Tree") %>%
                 # filter(taxon == "genus") %>%
                 # filter(cross_lab == "Within Lab") %>%
                 mutate(grouper = paste(classifier,predictee_lab,
                                        sep = "_",
                                        collapse = "_")),
               aes(x= predictor_lab,
                   y = misclass,
                   color = predictee_lab,
                   shape = cross_lab,
                   size = cross_lab,
                   alpha = cross_lab),
               size = .5) +
    geom_line(data = misclass %>%
                filter(transformation == trans) %>%
                filter(class_var == class_type) %>%
                # filter(taxon == "OTU") %>%
                # filter(classifier  == "Boosted Tree") %>%
                # filter(taxon %in% c("OTU","genus","order","phylum")) %>%
                filter(predictor_lab != predictee_lab) %>%
                mutate(grouper = paste(predictee_lab,
                                       sep = "_",
                                       collapse = "_")),
              aes(x = predictor_lab,
                  y = misclass,
                  color = predictee_lab,
                  group = predictee_lab,
                  alpha = cross_lab),
              size= .3
    ) +
    scale_color_brewer(palette = "Dark2") +
    scale_size_discrete(range = c(1,2)) +
    scale_alpha_discrete(range = c(.75,1)) +
    coord_cartesian(ylim= c(0,1)) +
    ylab("Misclassification") +
    xlab("Training Laboratory") +
    facet_grid(classifier ~ taxon) +
    # ggtitle("Proportion-Scale Specimen Classifier Performance by Taxon") +
    guides(color = guide_legend(title = "Test Laboratory", nrow = 2, byrow = 2, title.position = "top"),
           shape = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
                                override.aes = list(linetype  = 0)),
           size = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
                               override.aes = list(linetype  = 0)),
           alpha = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top",
                                override.aes = list(linetype  = 0))) +
    theme_bw(base_size = 8) +
    theme(axis.text.x=element_text(angle = 45, hjust = 1, size = 5),
          legend.position = "bottom")
  
  return(gg)
}
### Supplementary results figures

for(trans in c("proportion","proportion_diff","clr_pseudo_1","clr_pseudo_1_diff","presence")){
  for(cv in c("Specimen","Specimen Type")){
    
  fig_name <- paste("supp",trans,cv,".pdf",sep = "_",collapse = "_")
  
  pdf(width = 6, height = 3, file = fig_name)
    make_supp_results_fig(misclass,
                      trans,
                      cv) %>%
  print()
  dev.off()
  }
}


### (Supplementary) Discussion Overall Figure Including Uncentered Proportion and CLR
misclass_pretty <- misclass

misclass_pretty$transformation %<>% sapply(function(x) ifelse(x == "clr_pseudo_1_diff","CLR (Centered)",
                                                              ifelse(x == "clr_pseudo_1","CLR (Uncentered)",
                                                                     ifelse(x == "proportion_diff","Proportion (Centered)",
                                                                            ifelse(x == "proportion","Proportion (Uncentered)",
                                                                                   "Presence-Absence")))))

misclass_pretty$transformation  <- factor(misclass_pretty$transformation,
                                          levels = c("Proportion (Uncentered)",
                                                     "Proportion (Centered)",
                                                     "CLR (Uncentered)",
                                                     "CLR (Centered)",
                                                     "Presence-Absence"))



misclass_pretty$grouper <- 
  apply(misclass_pretty,1, function(x) paste(x["cross_lab"],x["classifier"], sep = "_", collapse = "_"))

pdf(width = 6, height = 3, file = "supp_median_misclass.pdf")
misclass_pretty %>%
  group_by %>%
  group_by(transformation, class_var, taxon, cross_lab, classifier) %>%
  summarize(median_misclass = median(misclass),
            grouper = unique(grouper))%>%
  
  ggplot() + 
  # geom_ribbon(aes(x = taxon, ymin = lower, ymax = upper, group = grouper, color = cross_lab,
  # fill = cross_lab), alpha = .5)+
  # geom_jitter(aes(x = taxon, y = misclass, color = cross_lab, shape = classifier), alpha = .2,
  #             height = 0, width = .3,
  # data = misclass_pretty) +
  geom_line(aes(x = taxon, y = misclass, color = cross_lab, linetype = classifier,
                group = interaction(predictee_lab, predictor_lab, classifier)), alpha = .2, size = .2,
            data = misclass_pretty) +
  geom_line(aes(x = taxon, y = median_misclass, group = interaction(cross_lab, classifier), color = cross_lab,
                linetype = classifier),
            alpha = 1,
            size = .5) + 
  facet_grid(class_var ~ transformation) + 
  theme_bw(base_size = 8) + 
  # ggtitle("Classification Results by Transformation, Taxon, and Classification Task") +
  ylab("Misclassification Error") + 
  xlab("Taxon") + 
  scale_shape_manual(values = c(4,1))+
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        text = element_text(size = 15),
        legend.position = "bottom"
  ) +
  guides(color = guide_legend(title = "Type of Prediction",nrow = 1, byrow = 1, title.position="top"),
         linetype = guide_legend(title = "Classifier",nrow = 1, byrow = 1, title.position="top"),
         fill = guide_legend(title = "Type of Prediction",nrow = 1, byrow = 1, title.position="top"),
         shape = guide_legend(title = "Classifier",nrow = 1, byrow = 1, title.position="top"))
dev.off()



misclass_pretty %<>%
  mutate(Centered = sapply(transformation, function(x) grepl("Centered",x,fixed = TRUE) %>% ifelse("Centered","Uncentered"))) %<>% 
  mutate(base_transformation = sapply(transformation, function(x) strsplit(as.character(x)," ", fixed = TRUE) %>% 
                                        unlist %>% (function(y) y[1])))
pdf(width = 4, height = 3, file = "supp_median_misclass_prop_clr.pdf")
misclass_pretty %>%
  filter(transformation != "Presence-Absence") %>%
  filter(class_var != "Specimen Type") %>%
  group_by(base_transformation, Centered, taxon, cross_lab, classifier) %>%
  summarize(median_misclass = median(misclass),
            q3 = quantile(misclass,.75),
            q1 = quantile(misclass,.25),
            grouper = unique(grouper))%>%
  ggplot() + 
  # geom_ribbon(aes(x = taxon, ymin = lower, ymax = upper, group = grouper, color = cross_lab,
  # fill = cross_lab), alpha = .5)+
  # geom_jitter(aes(x = taxon, y = misclass, color = cross_lab, shape = classifier), alpha = .2,
  #             height = 0, width = .3,
  # data = misclass_pretty) +
  # geom_line(aes(x = taxon, y = misclass, color = cross_lab, linetype = Centered,
  #               group = interaction(predictee_lab, predictor_lab, classifier)), alpha = .2, size = .2,
  #           data = misclass_pretty %>% filter(transformation != "Presence-Absence") %>%
  #             filter(class_var != "Specimen Type")) +
  geom_line(aes(x = taxon, y = median_misclass, group = interaction(cross_lab, Centered), color = cross_lab,
                linetype = Centered),
            alpha = 1,
            size = .1) + 
  
  geom_errorbar(aes(x = taxon, ymin = q1,ymax = q3, group = interaction(cross_lab, Centered), color = cross_lab,
                linetype = Centered),
            alpha = .6,
            width = .5,
            size = .1,
            position = position_dodge(.2)) + 
  facet_grid(base_transformation~classifier) + 
  theme_bw(base_size = 8) + 
  # ggtitle("Classification Results by Transformation, Taxon, and Classification Task") +
  ylab("Misclassification Error") + 
  xlab("Taxon") + 
  scale_linetype_manual(values = c(4,1))+
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
  ) +
  guides(color = guide_legend(title = "Type of Prediction",nrow = 1, byrow = 1, title.position="top"),
         linetype = guide_legend(title = "Centering",nrow = 1, byrow = 1, title.position="top"))
dev.off()




############################## Supplementary Tables ##############################


### make completeness tables

#calculate proportion of samples present (out of total run - should be multiple of 53 for each sequencing lab)
counts_non_preextracted_table <- metadata %>%
  filter(sequencing_wetlab %in% wet_labs) %>% #include only samples from valid sequencing labs
  filter(extraction_wetlab %in% extraction_wetlabs) %>% #include only samples from valid extraction labs
  filter(pre.extracted_bool == FALSE) %>% #include only non-pre-extracted samples
  # filter(as.character(sequencing_wetlab_archived) == #include only samples extracted and sequenced
  #          as.character(extraction_wetlab)) %>% #at same lab
  with(table(dry_lab,sequencing_wetlab))  #tabulate by dry lab and sequencing lab


pdf(width = 6,height= 2,"aliquot_counts.pdf")
counts_non_preextracted_table %>%
  (function(x) x[,!(colnames(x) %in% c("HL-F","HL-N"))]) %>%
  (function(x) grid.table(d = x, theme =  ttheme_minimal(base_size = 8,
                                                         padding = unit(c(1,2),"mm"))))
dev.off()

#determine which combinations of sequencing and dry labs
#(for samples extracted at sequencing lab only)
#have all specimens present
metadata$sequencing_wetlab %<>% as.factor()

pdf(width = 6,height= 2,"unique_specimens.pdf")
metadata %>%
  filter(sequencing_wetlab %in% wet_labs) %>% #include only samples from valid sequencing labs
  filter(extraction_wetlab %in% extraction_wetlabs) %>% #include only samples from valid extraction labs
  filter(pre.extracted_bool == FALSE) %>% #include only non-pre-extracted samples
  # filter(as.character(sequencing_wetlab) == #include only samples extracted and sequenced
  #          as.character(extraction_wetlab)) %>% #at same lab
  group_by(sequencing_wetlab, dry_lab, specimen) %>%
  summarize(unique_type = length(specimen_type_collapsed)) %>% #summarize some irrelevant variable
  with(table(dry_lab,sequencing_wetlab,specimen)) %>% #tabulate co-occurrences of dry lab, wet lab, and specimens
  apply(c(1,2), sum)  %>% #sum over specimen to get # unique specimens by dry lab wet lab combo
  (function(x) x[,!(colnames(x) %in% c("HL-F","HL-N"))]) %>%
  (function(x) grid.table(d = x, theme =  ttheme_minimal(base_size = 8,
                                                         padding = unit(c(1,2),"mm"))))
dev.off()


############################## Spikes in HL-B and HL-L misclassification for glmnet centered proportion data ##############################

# get confusion 
proportion_diff_conf <- get_confusion("glmnet",
              "proportion_diff")
proportion_diff_conf_info <- lapply(proportion_diff_conf, function(x)  data.frame(predictor_lab = x$predictor_lab,
                                                    class_var = x$class_var,
                                                    taxon = x$taxon,
                                                    classifier = x$classifier)) %>%
  (function(x) do.call(rbind,x))


### HLB Diagnostics

hlb_index <- which(with(proportion_diff_conf_info, predictor_lab == "HL-B" & class_var == "specimen" & taxon == "order"))

pdf(width = 13,height= 4,"confusion_hlb_order.pdf")
proportion_diff_conf[[hlb_index]]$confusion_matrices[["HL-B"]]  %>%
 (function(x) grid.table(d = x, theme =  ttheme_minimal(base_size = 8,
                                                        padding = unit(c(1,2),"mm"))))
dev.off()

hlb_order_test <- readRDS("HL-B_test_order_54332")

hlb_order_test[,sapply(colnames(hlb_order_test), function(x) grepl("k__",x,fixed = TRUE))] %<>%
  transform_data(transformation = "proportion_diff")

hlb_order_train <- readRDS("HL-B_train_order_54332")

hlb_order_train[,sapply(colnames(hlb_order_train), function(x) grepl("k__",x,fixed = TRUE))] %<>%
  transform_data(transformation = "proportion_diff") 

specimens <- hlb_order_train$specimen %>% unique() %>% (function(x) x[order(x)])

hlb_test_train_order_uncentered <- lapply(unique(hlb_order_test$specimen), function(z)
  rbind(
    hlb_order_test %>% filter(specimen == z) %>%
      select(starts_with("k__"), Bioinformatics.ID,dry_lab) %>%
      pivot_longer(starts_with("k__")) %>%
      mutate(group = "test",
             specimen = z),
    hlb_order_train %>% filter(specimen == z) %>%
      select(starts_with("k__"), Bioinformatics.ID,
             dry_lab) %>%
      pivot_longer(starts_with("k__")) %>%
      mutate(group = "train", 
             specimen = z))) %>%
  (function(x) do.call(rbind,x)) 


hlb_order_test[,sapply(colnames(hlb_order_test), function(x) grepl("k__",x,fixed = TRUE))] %<>%
  (function(x) center_data(x,specimens = hlb_order_test$specimen))


hlb_order_train[,sapply(colnames(hlb_order_train), function(x) grepl("k__",x,fixed = TRUE))] %<>%
  (function(x) center_data(x,specimens = hlb_order_train$specimen))

specimens <- hlb_order_train$specimen %>% unique() %>% (function(x) x[order(x)])

hlb_test_train_order <- lapply(unique(hlb_order_test$specimen), function(z)
  rbind(
    hlb_order_test %>% filter(specimen == z) %>%
      select(starts_with("k__"), Bioinformatics.ID,dry_lab) %>%
      pivot_longer(starts_with("k__")) %>%
      mutate(group = "test",
             specimen = z),
    hlb_order_train %>% filter(specimen == z) %>%
      select(starts_with("k__"), Bioinformatics.ID,
             dry_lab) %>%
      pivot_longer(starts_with("k__")) %>%
      mutate(group = "train", 
             specimen = z))) %>%
  (function(x) do.call(rbind,x)) 

hlb_glmnet_order <- readRDS(paste("registry_glmnet_train_proportion_diff/results/",hlb_index,".rds",sep = "",collapse = ""))

#pull lambda used for HL-B order classification
prop_diff_glmnet_params <- readRDS("glmnet_cv_results_proportion_diff")

hlb_order_lambda <-prop_diff_glmnet_params %>%
  filter(wet_lab == "HL-B", 
       taxon == "order",
       class_var == "specimen") %>%
  group_by(wet_lab) %>%
  summarize(lambda = lambda[which.min(cv)])

#use lambda to pull appropriate betas from glmnet object for HL-B order classification
hlb_order_betas <- hlb_glmnet_order$glmnet_results$beta %>% 
  lapply(function(x) 
    as.data.frame(as.matrix(x)[,
                 which.min( abs(as.numeric(hlb_glmnet_order$glmnet_results$lambda) - as.numeric(hlb_order_lambda)))]))

for(k in 1:length(hlb_order_betas)){
  hlb_order_betas[[k]]$specimen <- specimens[k]
  hlb_order_betas[[k]]$taxon <- rownames(hlb_order_betas[[k]])
}

hlb_order_betas %<>% (function(x) do.call(rbind,x))

colnames(hlb_order_betas)[1] <- "beta"
lin_preds <- vector(length(specimens), mode = "list")

for(k in 1:length(specimens)){
  print(k)
lp <- apply(hlb_test_train_order,1,
      (function(x) as.numeric(x["value"])*hlb_order_betas$beta[
      (hlb_order_betas$taxon == as.character(x["name"])) &
      (hlb_order_betas$specimen == specimens[k])]))
lin_preds[[k]] <- data.frame(lin_pred = lp,
                             specimen = specimens[k])
}

old_cols <- colnames(hlb_test_train_order)
for(k in 1:length(specimens)){
  hlb_test_train_order <- cbind(hlb_test_train_order, lin_preds[[k]]$lin_pred)
}

colnames(hlb_test_train_order) <- c(old_cols,specimens)


hlb_test_train_order %>%
  filter(dry_lab == "BL-2") %>%
  group_by(specimen) %>%
  summarize(length(unique(Bioinformatics.ID)))

lin_preds %<>% (function(x) do.call(rbind,x))
  
lps_hlb <- hlb_test_train_order %>%
  pivot_longer(cols = starts_with("sample"),
               names_to = "specimen_predicted",
               values_to = "lin_pred") %>%
  group_by(name) %>%
  mutate(max_abs_lp = max(abs(lin_pred))) %>%
  ungroup %>%
  filter(max_abs_lp >0)

pdf(width = 6, height = 6, "hlb_diagnostics.pdf")
lps_hlb %>%
  group_by(group,specimen,specimen_predicted,name) %>%
  summarize(lin_pred = mean(lin_pred)) %>%
  pivot_wider(names_from = group, values_from =lin_pred,
              id_cols = c(specimen,specimen_predicted,name)) %>%
  mutate("pred_63" = specimen_predicted == "sample63") %>%
  ggplot() + 
  geom_point(aes(x = test, y = train,
                 shape = pred_63,
                 color = sapply(name, function(x) ifelse(grepl("Lactobacillales",x),
                                                         "Lactobacillales",
                                                         ifelse(grepl("Sphingomonadales",x),
                                                                "Sphingomonadales",
                                                                "Other"))))) + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  facet_wrap(~specimen, scales = "free") +
  guides(shape = guide_legend(title = "Component of Linear \nPredictor for Sample 63",
                              title.position = "top"),
         color = guide_legend(title = "Order",
                              title.position = "top"))+
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 6)) +
  ylab("Linear Predictor Components on Training Set") +
  xlab("Linear Predictor Components on Test Set")
dev.off()

### HL-L Diagnostics

hll_index <- which(with(proportion_diff_conf_info, predictor_lab == "HL-L" & class_var == "specimen" & taxon == "genus"))

proportion_diff_conf[[hll_index]]$confusion_matrices[["HL-L"]] %>%
  apply(2, sum)

pdf(width = 13,height= 4,"confusion_hll_genus.pdf")
proportion_diff_conf[[hll_index]]$confusion_matrices[["HL-L"]]  %>%
  (function(x) grid.table(d = x, theme =  ttheme_minimal(base_size = 8,
                                                         padding = unit(c(1,2),"mm"))))
dev.off()


hll_genus_test <- readRDS("HL-L_test_genus_54332")

hll_genus_test[,sapply(colnames(hll_genus_test), function(x) grepl("k__",x,fixed = TRUE))] %<>%
  transform_data(transformation = "proportion_diff") 

hll_genus_train <- readRDS("HL-L_train_genus_54332")

hll_genus_train[,sapply(colnames(hll_genus_train), function(x) grepl("k__",x,fixed = TRUE))] %<>%
  transform_data(transformation = "proportion_diff") 

hll_test_train_genus_uncentered <- lapply(unique(hll_genus_test$specimen), function(z)
  rbind(
    hll_genus_test %>% filter(specimen == z) %>%
      select(starts_with("k__"), Bioinformatics.ID,dry_lab) %>%
      pivot_longer(starts_with("k__")) %>%
      mutate(group = "test",
             specimen = z),
    hll_genus_train %>% filter(specimen == z) %>%
      select(starts_with("k__"), Bioinformatics.ID,
             dry_lab) %>%
      pivot_longer(starts_with("k__")) %>%
      mutate(group = "train", 
             specimen = z))) %>%
  (function(x) do.call(rbind,x)) 


hll_genus_test[,sapply(colnames(hll_genus_test), function(x) grepl("k__",x,fixed = TRUE))] %<>%
  (function(x) center_data(x,specimens = hll_genus_test$specimen))


hll_genus_train[,sapply(colnames(hll_genus_test), function(x) grepl("k__",x,fixed = TRUE))] %<>%
  (function(x) center_data(x,specimens = hll_genus_train$specimen))


specimens <- hll_genus_train$specimen %>% unique() %>% (function(x) x[order(x)])

hll_test_train_genus <- lapply(unique(hll_genus_test$specimen), function(z)
  rbind(
    hll_genus_test %>% filter(specimen == z) %>%
      select(starts_with("k__"), Bioinformatics.ID,dry_lab) %>%
      pivot_longer(starts_with("k__")) %>%
      mutate(group = "test",
             specimen = z),
    hll_genus_train %>% filter(specimen == z) %>%
      select(starts_with("k__"), Bioinformatics.ID,
             dry_lab) %>%
      pivot_longer(starts_with("k__")) %>%
      mutate(group = "train", 
             specimen = z))) %>%
  (function(x) do.call(rbind,x)) 



hll_glmnet_genus <- readRDS(paste("registry_glmnet_train_proportion_diff/results/",hll_index,".rds",sep = "",collapse = ""))

#pull lambda used for HL-B order classification
prop_diff_glmnet_params <- readRDS("glmnet_cv_results_proportion_diff")

hll_genus_lambda <-prop_diff_glmnet_params %>%
  filter(wet_lab == "HL-L", 
         taxon == "genus",
         class_var == "specimen") %>%
  group_by(wet_lab) %>%
  summarize(lambda = lambda[which.min(cv)])

#use lambda to pull appropriate betas from glmnet object for HL-B order classification
hll_genus_betas <- hll_glmnet_genus$glmnet_results$beta %>% 
  lapply(function(x) 
    as.data.frame(as.matrix(x)[,
                               which.min( abs(as.numeric(hll_glmnet_genus$glmnet_results$lambda) - as.numeric(hll_genus_lambda)))]))

for(k in 1:length(hll_genus_betas)){
  hll_genus_betas[[k]]$specimen <- specimens[k]
  hll_genus_betas[[k]]$taxon <- rownames(hll_genus_betas[[k]])
}

hll_genus_betas %<>% (function(x) do.call(rbind,x))

colnames(hll_genus_betas)[1] <- "beta"
lin_preds <- vector(length(specimens), mode = "list")

for(k in 1:length(specimens)){
  print(k)
  lp <- apply(hll_test_train_genus,1,
              (function(x) as.numeric(x["value"])*hll_genus_betas$beta[
                (hll_genus_betas$taxon == as.character(x["name"])) &
                  (hll_genus_betas$specimen == specimens[k])]))
  lin_preds[[k]] <- data.frame(lin_pred = lp,
                               specimen = specimens[k])
}


old_cols <- colnames(hll_test_train_genus)
for(k in 1:length(specimens)){
  hll_test_train_genus <- cbind(hll_test_train_genus, lin_preds[[k]]$lin_pred)
}

colnames(hll_test_train_genus) <- c(old_cols,specimens)


lps_hll <- hll_test_train_genus %>%
  pivot_longer(cols = starts_with("sample"),
               names_to = "specimen_predicted",
               values_to = "lin_pred") %>%
  group_by(name) %>%
  mutate(max_abs_lp = max(abs(lin_pred))) %>%
  ungroup %>%
  filter(max_abs_lp >0)

pdf(width = 6, height = 6, "hll_diagnostics.pdf")
lps_hll %>%
  group_by(group,specimen,specimen_predicted,name) %>%
  summarize(lin_pred = mean(lin_pred)) %>%
  pivot_wider(names_from = group, values_from =lin_pred,
              id_cols = c(specimen,specimen_predicted,name)) %>%
  mutate("pred_ofint" = sapply(specimen_predicted, 
                               function(x) ifelse(x == "sample101","Sample 101", 
                                                  ifelse(x == "sample61","Sample 61", "Other Sample")))) %>%
  ggplot() + 
  geom_point(aes(x = test, y = train,
                 color = pred_ofint)) + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  facet_wrap(~specimen, scales = "free") +
  guides(color = guide_legend(title = "Component of Linear \nPredictor for Which Sample?"))+
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 6)) +
  ylab("Linear Predictor Components on Training Set") +
  xlab("Linear Predictor Components on Test Set")
dev.off()



######## Sensitivity Analysis: Classification using Dry Labs Individually ########

### get sensitivity analysis results
dl_results <- vector(8, mode = "list")

### sensitivity analyses using data from dry labs individually
included_wet_labs_archive <- included_wet_labs
#exclude wet labs with only one aliquot set for sensitivity analysis
included_wet_labs <- included_wet_labs[!(included_wet_labs %in% 
                                           c("HL-I","HL-K","HL-L"))]
glmnet_params_archive <- glmnet_params
glmnet_params <- glmnet_params %>% dplyr::filter(!(Var1 %in% 
                                                     c("HL-I","HL-K","HL-L"))) %>%
  filter(Var3 == "specimen")

xgb_params_archive <- xgb_params
xgb_params <- xgb_params %>% dplyr::filter(!(Var1 %in% 
                                               c("HL-I","HL-K","HL-L"))) %>%
  filter(Var6 == "specimen")

counter <- 0
for(dl in included_dry_labs){
  for(transformation in c("proportion_diff","clr_pseudo_1_diff")){
    counter <- counter + 1
    temp_result <- get_misclass_results("glmnet",
                                        paste(transformation,dl,sep = "_",
                                              collapse = "_"),
                                        nresults = 70)
    temp_result$dry_lab <- dl
    dl_results[[counter]] <- temp_result
  }
}

dl_results <- do.call(rbind,dl_results)
dl_results$transformation <- sapply(dl_results$transformation,
                                    function(x) ifelse(grepl("proportion",
                                                             x,
                                                             fixed = TRUE),
                                    "Proportion (Centered)",
                                    "CLR (Centered)"))

dl_results$same_lab <- as.character(dl_results$predictor_lab) ==
  as.character(dl_results$predictee_lab)
dl_results$same_lab <- sapply(dl_results$same_lab,
                              function(x) ifelse(x,
                                                 "Within Laboratory",
                                                 "Across Laboratory"))

dl_results$traintest <- apply(dl_results,1,
                              function(x) paste(c(x["predictor_lab"],
                                                x["predictee_lab"]),
                                                sep = "_",
                                                collapse = "_"))

dl_results$predictor_lab %<>% factor(levels = c("HL-J",
                                                "HL-H",
                                                "HL-F_1",
                                                "HL-B",
                                                "HL-C"))
levels(dl_results$predictor_lab)[levels(dl_results$predictor_lab) == "HL-F_1"] <- 
  "HL-F"

misclass_dl <- misclass %>%
  filter(classifier == "Elastic Net") %>%
  filter(predictee_lab %in% c("HL-J","HL-H", "HL-F","HL-B",
                              "HL-C")) %>%
  filter(predictor_lab %in% c("HL-J","HL-H", "HL-F","HL-B",
                              "HL-C")) %>%
  filter(class_var == "Specimen") %>%
  filter(transformation %in% c("proportion_diff",
                               "clr_pseudo_1_diff")) %>%
  mutate(dry_lab = "All") %>%
  mutate(same_lab = sapply(cross_lab,
                           function(x) ifelse(x == "Within Lab",
                                              "Within Laboratory",
                                              "Across Laboratory"))) %>%
  select(-cross_lab)

misclass_dl$class_var <- "specimen"

misclass_dl$transformation <- sapply(misclass_dl$transformation,
                                    function(x) ifelse(grepl("proportion_diff",
                                                             x,
                                                             fixed = TRUE),
                                                       "Proportion (Centered)",
                                                       "CLR (Centered)"))

misclass_dl <- rbind(misclass_dl %>% mutate(dry_lab = "BL-2"),
                     misclass_dl %>% mutate(dry_lab = "BL-4"),
                     misclass_dl %>% mutate(dry_lab = "BL-6"),
                     misclass_dl %>% mutate(dry_lab = "BL-8")
                     )

misclass_dl$traintest <- apply(misclass_dl,1,
                              function(x) paste(c(x["predictor_lab"],
                                                  x["predictee_lab"]),
                                                sep = "_",
                                                collapse = "_"))

misclass_dl$predictor_lab %<>% factor(levels = c("HL-J",
                                                "HL-H",
                                                "HL-F",
                                                "HL-B",
                                                "HL-C"))

misclass_dl$predictee_lab %<>% factor(levels = c("HL-J",
                                                 "HL-H",
                                                 "HL-F",
                                                 "HL-B",
                                                 "HL-C"))

dl_results$type <- "Single Bioinformatics"
misclass_dl$type <- "All Bioinformatics"
dl_results <- rbind(dl_results, misclass_dl)

dl_results$predictee_lab[dl_results$predictee_lab == "HL-F_1"] <- "HL-F"

dl_results$type %<>% factor(levels = c(
                            "All Bioinformatics","Single Bioinformatics"))

##what's going on?

pdf(width = 6, height = 4.5, file = "dry_proportion.pdf")
ggplot() +
  geom_line(data = dl_results %>%
              filter(transformation %in% c("Proportion (Centered)")) %>%
              filter(class_var == "specimen") %>%
              # filter(taxon == "OTU") %>%
              # filter(classifier  == "Boosted Tree") %>%
              filter(taxon %in% c("OTU","genus","order","phylum")) %>%
              filter(predictor_lab != predictee_lab),
            aes(x = predictor_lab,
                y = misclass,
                color = predictee_lab,
                group = interaction(predictee_lab,type),
                alpha = type,
                linetype = type),
            size= .3) +
  geom_point(data = dl_results %>%
               filter(transformation %in% c("Proportion (Centered)")) %>%
               filter(taxon %in% c("OTU","genus","order","phylum")) %>%
               filter(class_var == "specimen") %>%
               # filter(classifier == "Boosted Tree") %>%
               # filter(taxon == "genus") %>%
               # filter(cross_lab == "Within Lab") %>%
               mutate(grouper = paste(classifier,predictee_lab,
                                      sep = "_",
                                      collapse = "_")),
             aes(x= predictor_lab,
                 y = misclass,
                 color = predictee_lab,
                 shape = same_lab,
                 size = same_lab,
                 alpha = type),
             size = .5) +
  scale_color_brewer(palette = "Dark2") +
  scale_size_discrete(range = c(1,2)) +
  scale_alpha_discrete(range = c(.5,1)) +
  scale_linetype_manual(values = c(2,1)) +
  coord_cartesian(ylim= c(0,1)) +
  ylab("Misclassification") +
  xlab("Training Laboratory") +
  facet_grid(dry_lab ~ taxon) +
  # ggtitle("Proportion-Scale Specimen Classifier Performance by Taxon") +
  guides(color = guide_legend(title = "Test Laboratory", nrow = 4, byrow = 4, title.position = "top"),
         shape = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top"),
         size = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top"),
         alpha = guide_legend(title = "Training and Validation", nrow = 2, byrow = 2, title.position = "top"),
         linetype =guide_legend(title = "Training and Validation", nrow = 2, byrow = 2, title.position = "top") ) +
  theme_bw(base_size = 8) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 8))
dev.off()

pdf(width = 6, height = 4.5, file = "dry_clr.pdf")
ggplot() +
  geom_line(data = dl_results %>%
              filter(transformation %in% c("CLR (Centered)")) %>%
              filter(class_var == "specimen") %>%
              # filter(taxon == "OTU") %>%
              # filter(classifier  == "Boosted Tree") %>%
              filter(taxon %in% c("OTU","genus","order","phylum")) %>%
              filter(predictor_lab != predictee_lab),
            aes(x = predictor_lab,
                y = misclass,
                color = predictee_lab,
                group = interaction(predictee_lab,type),
                alpha = type,
                linetype = type),
            size= .3) +
  geom_point(data = dl_results %>%
               filter(transformation %in% c("CLR (Centered)")) %>%
               filter(taxon %in% c("OTU","genus","order","phylum")) %>%
               filter(class_var == "specimen") %>%
               # filter(classifier == "Boosted Tree") %>%
               # filter(taxon == "genus") %>%
               # filter(cross_lab == "Within Lab") %>%
               mutate(grouper = paste(classifier,predictee_lab,
                                      sep = "_",
                                      collapse = "_")),
             aes(x= predictor_lab,
                 y = misclass,
                 color = predictee_lab,
                 shape = same_lab,
                 size = same_lab,
                 alpha = type),
             size = .5) +
  scale_color_brewer(palette = "Dark2") +
  scale_size_discrete(range = c(1,2)) +
  scale_alpha_discrete(range = c(.5,1)) +
  scale_linetype_manual(values = c(2,1)) +
  coord_cartesian(ylim= c(0,1)) +
  ylab("Misclassification") +
  xlab("Training Laboratory") +
  facet_grid(dry_lab ~ taxon) +
  # ggtitle("Proportion-Scale Specimen Classifier Performance by Taxon") +
  guides(color = guide_legend(title = "Test Laboratory", nrow = 4, byrow = 4, title.position = "top"),
         shape = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top"),
         size = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top"),
         alpha = guide_legend(title = "Training and Validation", nrow = 2, byrow = 2, title.position = "top"),
         linetype =guide_legend(title = "Training and Validation", nrow = 2, byrow = 2, title.position = "top") ) +
  theme_bw(base_size = 8) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 8))
dev.off()


#### scaling sensitivity analysis and phylum-level CRC case-control comparison

included_wet_labs <- included_wet_labs_archive
glmnet_params <- glmnet_params_archive

misclass_glmnet <- lapply(c("proportion_diff", "proportion_diff_scale",
                            "clr_pseudo_1_diff", "clr_pseudo_1_diff_scale",
                            "presence"),
                          function(x) get_misclass_results("glmnet",x))

misclass_glmnet %<>% (function(x) do.call(rbind,x))


misclass <- rbind(misclass_glmnet)

misclass$cross_lab <- (as.character(misclass$predictor_lab) ==
                         as.character(misclass$predictee_lab)) %>%
  sapply(function(x) ifelse(x,"Within Lab","Cross Lab"))
# #
misclass$class_var[misclass$class_var == "specimen"] <- "Specimen"
misclass$class_var[misclass$class_var == "specimen_type_collapsed"] <- "Specimen Type"
levels(misclass$class_var) <- c("Specimen", "Specimen Type")

misclass$predictee_lab %<>% as.factor()
misclass$predictor_lab %<>% as.factor()

levels(misclass$predictor_lab)[levels(misclass$predictor_lab)  == "HL-F_1"] <- "HL-F"
levels(misclass$predictee_lab)[levels(misclass$predictee_lab)  == "HL-F_1"] <- "HL-F"

scale.names <- c('Elastic Net.proportion_diff' = "Proportion \nCentered",
                   'Elastic Net.proportion_diff_scale' = "Proportion \nCentered, Scaled",
                   'Elastic Net.clr_pseudo_1_diff' = "CLR \nCentered",
                   'Elastic Net.clr_pseudo_1_diff_scale' = "CLR \nCentered, Scaled",
                   'OTU' = "OTU",
                   'genus' = "Genus",
                   'order' = "Order",
                   'phylum' = "Phylum")

misclass$trans <- sapply(misclass$transformation,
                         function(x) ifelse(grepl("proportion",x),"Proportion (Centered)",
                                                  "CLR (Centered)"))

misclass$scaled <- sapply(misclass$transformation,
                          function(x) ifelse(grepl("scale",x),"Rescaled",
                                             "Not Rescaled"))

pdf(width = 6, height = 4, file = "scale_comparison.pdf")
ggplot() +
  geom_line(data = misclass %>%
              filter(transformation %in% c("proportion_diff","clr_pseudo_1_diff",
                                           "proportion_diff_scale",
                                           "clr_pseudo_1_diff_scale")) %>%
              filter(class_var == "Specimen") %>%
              # filter(taxon == "OTU") %>%
              # filter(classifier  == "Boosted Tree") %>%
              filter(taxon %in% c("OTU","genus","order","phylum")) %>%
              filter(predictor_lab != predictee_lab) %>%
              mutate(grouper = paste(predictee_lab,
                                     sep = "_",
                                     collapse = "_")),
            aes(x = predictor_lab,
                y = misclass,
                color = predictee_lab,
                group = interaction(predictee_lab, scaled)),
            size= .3
  ) +
  geom_point(data = misclass %>%
               filter(transformation %in% c("proportion_diff","clr_pseudo_1_diff",
                                            "proportion_diff_scale",
                                            "clr_pseudo_1_diff_scale")) %>%
               filter(class_var == "Specimen") %>%
               filter(taxon %in% c("OTU","genus","order","phylum")) %>%
               # filter(classifier == "Boosted Tree") %>%
               # filter(taxon == "genus") %>%
               # filter(cross_lab == "Within Lab") %>%
               mutate(grouper = paste(classifier,predictee_lab,
                                      sep = "_",
                                      collapse = "_")),
             aes(x= predictor_lab,
                 y = misclass,
                 color = predictee_lab,
                 shape = cross_lab,
                 size = cross_lab),
             size = .5) +
  scale_color_brewer(palette = "Dark2") +
  scale_size_discrete(range = c(1,2)) +
  scale_alpha_manual(values = c(.5,1)) +
  scale_linetype_manual(values = c(2,1))+
  coord_cartesian(ylim= c(0,1)) +
  ylab("Misclassification") +
  xlab("Training Laboratory") +
  facet_grid(interaction(classifier,transformation) ~ taxon,  
             labeller = as_labeller(scale.names)) +
  # ggtitle("Proportion-Scale Specimen Classifier Performance by Taxon") +
  guides(color = guide_legend(title = "Test Laboratory", nrow = 4, byrow = 4, title.position = "top"),
         shape = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top"),
         size = guide_legend(title = "Type of Prediction", nrow = 2, byrow = 2, title.position = "top")) +
  theme_bw(base_size = 8) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 8))
dev.off()

### CRC case/control comparisons for 4 phyla considered by Sinha et al.
############################## Phylum-level Plots (Supplementary) ############################## 
phylum_data <- readRDS("HL-B_train_phylum_54332")

phylum_data$transformation <- "Proportion"

phylum_data <- phylum_data[-(1:nrow(phylum_data)),]

for(wl in included_wet_labs){
  new_data <- rbind( readRDS(paste(wl,"_train_phylum_54332",
                                   sep = "", 
                                   collapse = "")),
                     readRDS(paste(wl,"_test_phylum_54332",
                                   sep = "",
                                   collapse = "")
                     ))
  new_clr_data <- new_data
  new_clr_data[,sapply(colnames(new_data),
                       function(x) grepl("k__",x,fixed = TRUE))] %<>%
    transform_data("clr_pseudo_1")
  new_clr_data$transformation <- "CLR"
  
  new_proportion_data <- new_data
  new_proportion_data[,sapply(colnames(new_data),
                              function(x) grepl("k__",x,fixed = TRUE))] %<>%
    transform_data("proportion") 
  new_proportion_data$transformation <- "Proportion"
  
  
  phylum_data <- rbind(phylum_data,
                       new_clr_data,
                       new_proportion_data)
}

phylum_data_long <- phylum_data %>%
  pivot_longer(starts_with("k__"))



phylum_data_long$specimen_type_collapsed %<>%
  factor(levels = c("Oral artificial colony",
                    "Fecal artificial colony",
                    "Freeze-dried",
                    "Fresh",
                    "Robogut"   
  ))
phylum_data_long$specimen %<>%
  factor(levels = c("sample10", #fresh
                    "sample101", #robogut
                    "sample55", #freeze-dried
                    "sample11", #fresh
                    "sample56", #freeze-dried
                    "sample12", #fresh
                    "sample57", #freeze-dried
                    "sample13", #fresh
                    "sample103", #oral mock
                    "sample104", #fecal mock
                    "sample14", #fresh
                    "sample59", #freeze-dried
                    "sample2", #fresh
                    "sample61", #freeze-dried
                    "sample3", #fresh
                    "sample63", #freeze-dried
                    "sample4", #fresh
                    "sample64", #freeze-dried
                    "sample6", #fresh
                    "sample8", #fresh,
                    "sample9", #fresh
                    "sample102" #robogut
  ))



mononym_conversion <- cbind(c("Firmicutes",
                              "Bacteroidetes",
                              "Proteobacteria",
                              "Actinobacteria"),
                            c("k__Bacteria|p__Firmicutes|",
                              "k__Bacteria|p__Bacteroidetes|",
                              "k__Bacteria|p__Proteobacteria|",
                              "k__Bacteria|p__Actinobacteria|")
)

phylum_data_long$health_status[is.na(phylum_data_long$health_status)] <- 
  as.character(phylum_data_long$specimen_type_collapsed)[is.na(phylum_data_long$health_status)]

phylum_data_long$health_status[phylum_data_long$specimen_type_collapsed == "Robogut"] <- "Robogut" 

phylum_data_long$health_status %<>% factor(
  levels = c("Control",
             "CRC case",
             "Sick (Diabetes/RA)",
             "Sick (ICU)",
             "Healthy",
             "Robogut",
             "Fecal artificial colony",
             "Oral artificial colony")
)

for(i in 1:4){
  
  pdf(paste("fresh_",mononym_conversion[i,1],".pdf",
            sep = "",collapse = ""),
      width = 6.5,
      height = 7)

  temp_plot <-
    phylum_data_long %>%
    filter(name %in% mononym_conversion[i,2]) %>% 
    filter(specimen_type_collapsed %in% c("Fresh")) %>% 
    ggplot() + 
    geom_jitter(aes(x = interaction(specimen,health_status), y = value, 
                    color = specimen,
                    group = interaction(specimen,
                                        Bioinformatics.ID),
                    shape = dry_lab),
                
                height =0,
                width = .125,
                size = .5) + 
    scale_color_viridis_d() +  
    facet_grid(transformation~sequencing_wetlab,
               # interaction(Summary,transformation)~sequencing_wetlab,
               scales = "free_y"#,
               # labeller = as_labeller(transformation_names)
    ) + 
    theme_bw() + 
    xlab("Specimen") + 
    ylab("Relative Abundance (CLR or Proportion Scale)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
          legend.position = "bottom",
          legend.box = "horizontal",
          text = element_text(size = 10),
          axis.text = element_text(size = 5),
          axis.title = element_text(size = 10),
          plot.title = element_text(size=10),
          panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")
    ) + 
    guides(
      color = guide_legend(title = "Specimen",
                           legend.title=element_text(size=8),
                           title.position="top",
                           ncol = 4),
      shape = guide_legend(title = "Bioinformatics \nLaboratory",
                           legend.title=element_text(size=8),
                           title.position="top",
                           ncol = 4)) +
    ggtitle(paste("Measured Abundance of", mononym_conversion[i,1], "among Fresh Human Specimens
by Specimen, Health status, Transformation, Sequencing Laboratory, \nand Bioinformatics Laboratory",
                  sep = " ",collapse =" "))
  print(temp_plot)
  dev.off()
  
}


