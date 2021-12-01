
paper_extractor <- function(dummy_var){ #function to extract mbqc paper data

  print(dummy_var)
  gc()
  if(file.exists("nbt.3981-S9.zip")){
    print("trying to read in data")

  } else{
            stop("File nbt.3981-S9.zip not found. \n
            To proceed, download supplementary data set 6 \n
            at https://www.nature.com/articles/nbt.3981 \n
            and rename downloaded file nbt.3981-S9.zip as necessary.")
                              }


  ncol_paper <- 28435


  nchunks <- ceiling(ncol_paper/5000)

  dir.create("paper_temp")

  for(chunk in 1:nchunks){
    # system('free -m')
    print(paste("reading in chunk ", chunk, sep = "", collapse = " "))

    temp <- read.table(file = paste(getwd(),"/nbt.3981-S9.zip",sep = ""),
                       header = FALSE,
                       sep = "\t",
                       skip =(chunk - 1)*5000,
                       nrows = 5000)
    cnames <- temp[,1]

    temp <- as.data.table(t(temp[,-1]))
    
    colnames(temp) <- cnames

    print("data loaded")

    saveRDS(temp,
            file = paste("paper_temp/paper_chunk",chunk, sep = "", collapse = "")
            )
  }


  mbqc_paper_data <- readRDS(file = paste("paper_temp/paper_chunk",1, sep = "",
                                                collapse = ""))


  # saveRDS(mbqc_paper_data, file = "mbqc_paper_data")
  # rm("mbqc_paper_data")

  for(chunk in 2:nchunks){
    print("reading in new chunk")
    temp <- readRDS(file = paste("paper_temp/paper_chunk",chunk, sep = "",
                                 collapse = ""))

    print("reading in mbqc_paper_data")
    # mbqc_paper_data <- readRDS("mbqc_paper_data")

    mbqc_paper_data <- cbind(mbqc_paper_data,temp)

    rm(temp)

    
    
    # if(chunk < nchunks){
    # rm(mbqc_paper_data)
    # }
    
    gc()

    print(paste("added chunk", chunk, collapse = ""))

  }
  
  


  

  unlink("paper_temp", recursive = TRUE)

  print("temporary data removed")
  
 

  # # properly format column names (currently first row)
  # colnames(mbqc_paper_data) <- mbqc_paper_data[1,]
  # mbqc_paper_data %<>% (function(x) x[-1,])

  # # format as data table
  # mbqc_paper_data %<>% as.data.table

  # format count data as numeric variables
  to_numeric <- sapply(colnames(mbqc_paper_data),
                       function(x) grepl("k__",x,fixed = TRUE))
  
  to_numeric <- colnames(mbqc_paper_data)[to_numeric]
  
  
  mbqc_paper_data[,(to_numeric):= lapply(.SD, as.numeric), .SDcols = to_numeric]
  mbqc_paper_data %<>% as.data.frame()
  


  # filter out rows with invalid values of sequencing wet lab
  mbqc_paper_data %<>%
    filter(!is.na(sequencing_wetlab))
  mbqc_paper_data %<>%
    filter(sequencing_wetlab != "Unknown")


  # filter out rows with invalid values of bioinformatics lab
  mbqc_paper_data %<>%
    filter(!is.na(dry_lab))
  mbqc_paper_data %<>%
    filter(dry_lab != "Unknown")
  mbqc_paper_data %<>%
    filter(dry_lab != "BL-10") #BL-10 not included in MBQC analysis

  # filter out rows with invalid values of specimen
  mbqc_paper_data %<>%
    filter(!is.na(specimen)) %<>%
    filter(specimen != "Unknown")

  # save full dataset
  if(!dir.exists("mbqc_paper_data_products")){
    dir.create("mbqc_paper_data_products")
  }
  saveRDS(mbqc_paper_data,
          file = "mbqc_paper_data_products/full_data_table")

  # save data sets for each sequencing lab
  seq_labs <- mbqc_paper_data$sequencing_wetlab %>% unique

  for(lab in seq_labs){
    mbqc_paper_data %>%
      dplyr::filter(sequencing_wetlab == lab) %>%
      saveRDS(file = paste("mbqc_paper_data_products/",lab,sep = ""))
    print(paste("saved data for wet lab ", lab,sep = "", collapse = ""))
  }

  # save metadata only
  mbqc_paper_data %>%
    dplyr::select(-starts_with("k__")) %>%
    saveRDS(file = "mbqc_paper_data_products/meta_only")
  
  print("saved metadata")



  return("done")}

### hmp site data extraction function
site_extractor <- function(dummy_var){
  
  print(dummy_var)
  gc()
  if(file.exists("mbqc_integrated_otus.tsv.gz")){
    print("trying to read in data")
 
  } else{
                                  stop("File mbqc_integrated_otus.tsv.gz not found. \n 
            To proceed, download integrated OTU table  \n
            at http://downloads.ihmpdcc.org/data/MBQC/mbqc_integrated_otus.tsv.gz \n
            and rename downloaded file mbqc_integrated_otus.tsv.gz as necessary.")
  }
  print("successfully read in data")
  # transpose object so variables are columns
  
  ncol_site <- 27282
  
  
  nchunks <- ceiling(ncol_site/5000)
  
  dir.create("site_temp")
  
  for(chunk in 1:nchunks){
    # system('free -m')
    print(paste("reading in chunk ", chunk, sep = "", collapse = " "))
    
    temp <- read.table(file = paste(getwd(),"/mbqc_integrated_otus.tsv.gz",sep = ""),
                       header = FALSE,
                       sep = "\t",
                       skip =(chunk - 1)*5000,
                       nrows = 5000)
    
    cnames <- temp[,1]
    
    temp <- as.data.table(t(temp[,-1]))
    
    colnames(temp) <- cnames
    
    print("data loaded")
    
    saveRDS(temp,
            file = paste("site_temp/site_chunk",chunk, sep = "", collapse = "")
    )
  }
  
  
  mbqc_site_data <- readRDS(file = paste("site_temp/site_chunk",1, sep = "",
                                          collapse = ""))
  
  # saveRDS(mbqc_paper_data, file = "mbqc_paper_data")
  # rm("mbqc_paper_data")
  
  for(chunk in 2:nchunks){
    print("reading in new chunk")
    temp <- readRDS(file = paste("site_temp/site_chunk",chunk, sep = "",
                                 collapse = ""))
    
    print("reading in mbqc_site_data")
    
    mbqc_site_data <- cbind(mbqc_site_data,temp)
    
    rm(temp)
    

    # }
    
    gc()
    
    print(paste("added chunk", chunk, collapse = ""))
    
  }
  
  
  
  
  print("data transposed")
  
  unlink("site_temp", recursive = TRUE)
  
  
  # properly format column names (currently first row)
  # colnames(mbqc_site_data) <- mbqc_site_data[1,]
  # mbqc_site_data %<>% (function(x) x[-1,])
  # 
  # # format as data table
  # mbqc_site_data %<>% as.data.table
  
  # format count data as numeric variables
  to_numeric <- sapply(colnames(mbqc_site_data),
                       function(x) grepl("k__",x,fixed = TRUE))
  
  to_numeric <- colnames(mbqc_site_data)[to_numeric]
  
  
  mbqc_site_data[,(to_numeric):= lapply(.SD, as.numeric), .SDcols = to_numeric]
  
  mbqc_site_data %<>% as.data.frame()
  
  # filter out rows with invalid values of sequencing wet lab
  mbqc_site_data %<>%
    filter(!is.na(sequencing_wetlab))
  mbqc_site_data %<>%
    filter(sequencing_wetlab != "Unknown")
  
  
  # filter out rows with invalid values of bioinformatics lab
  mbqc_site_data %<>%
    filter(!is.na(dry_lab))
  mbqc_site_data %<>%
    filter(dry_lab != "Unknown")
  mbqc_site_data %<>%
    filter(dry_lab != "BL-10") #BL-10 not included in MBQC analysis
  
  
  # save full dataset
  if(!dir.exists("hmp_mbqc_site_data_products")){
    dir.create("hmp_mbqc_site_data_products")
  }
  saveRDS(mbqc_site_data,
          file = "hmp_mbqc_site_data_products/full_data_table")
  
  # save data sets for each sequencing lab
  seq_labs <- mbqc_site_data$sequencing_wetlab %>% unique
  
  for(lab in seq_labs){
    mbqc_site_data %>%
      dplyr::filter(sequencing_wetlab == lab) %>%
      saveRDS(file = paste("hmp_mbqc_site_data_products/",lab,sep = ""))
    print(paste("saved data for wet lab ", lab,sep = "", collapse = ""))
  }
  
  # save metadata only
  
  mbqc_site_data %>%
    dplyr::select(-starts_with("k__")) %>%
    saveRDS(file = "hmp_mbqc_site_data_products/meta_only")
  
  print("saved metadata")
  
  return("done")
}

### taxonomic aggregation function
aggregate_level <- function(wet_lab_data, level){
  if(!(level %in% c("OTU","species","genus","family","order","class","phylum"))){
    stop("Invalid taxonomic level")
  }
  if(level == "OTU"){ # do nothing if we're looking at OTU level data
    return(wet_lab_data)
  }
  # get highest level to delete
  deletion_level <- c(NA,
             "OTU",
             "s__",
             "g__",
             "f__",
             "o__",
             "c__")[
               which(c("OTU","species","genus","family","order","class","phylum") == 
                       level)]
  ### get all taxa
  if(deletion_level == "OTU"){
  taxa <- wet_lab_data %>%
    select(starts_with("k__")) %>%
    colnames() %>%
      lapply(function(x) list(gregexpr("\\|",x)[[1]] %>% last,x)) %>%
      sapply(function(x) substr(x[[2]],1,x[[1]]-1))
  } else{
  
  taxa <- wet_lab_data %>%
    select(starts_with("k__")) %>%
    colnames() %>% 
    lapply(function(x) 
      list(gregexpr(deletion_level,x)[[1]][1],
                                  x)) %>%
    sapply(function(x) substr(x[[2]],start = 1, stop = x[[1]]-1))
  }
  
  # get unique taxa
  
  unique_taxa <- taxa %>% unique
  
  
  # get aggregated data
  
  aggregated_data <- data.frame(matrix(nrow = nrow(wet_lab_data),
                    ncol = length(unique_taxa)))
  
  aggregated_data %<>% apply(2, as.numeric)
  
  colnames(aggregated_data) <- unique_taxa
  
  #save metadata
  metadata <- wet_lab_data %>%
    select(-starts_with("k__"))
  
  #separate count data
  wet_lab_data %<>%
    select(starts_with("k__"))
  
  #aggregate
  n_taxa <- length(unique_taxa)
  for(i in 1:length(unique_taxa)){
    print(paste(level, i, "of", n_taxa))
    aggregated_data[,i] <- 
        wet_lab_data %>%
        select(which(taxa == unique_taxa[i])) %>%
        apply(2, as.numeric) %>%
        apply(1, sum)
  }
  
  aggregated_data <- cbind(metadata,aggregated_data)
  
  return(aggregated_data)
}

### function to split out test and training sets at a given level of taxonomic aggregation
split_test <- function(wet_lab,
                       included_dry_labs,
                       seed = 54332,
                       level){

  
  #print wet lab
  print(wet_lab)
  wet_lab %<>% as.character
  wl_name_length <- wet_lab %>% nchar
  
  if(substr(wet_lab,wl_name_length,wl_name_length) %in% c("1","2")){
    wet_sublab <- substr(wet_lab,wl_name_length, wl_name_length)
    wet_lab <- substr(wet_lab,1,(wl_name_length -2))
    sublab_ind <- TRUE
  } else{
    sublab_ind <- FALSE
    wet_sublab <- ""
  }
  
  #convert to name in wet lab files
  wet_lab_file <- sub("HL-",x=wet_lab,replacement="")
  wet_lab_file <- paste("HL-",wet_lab_file,sep ="")
  
  #load all data for wet lab
  # load(paste("wet_lab_data_files/wet_lab_",wet_lab_file,".RData", sep = ""))
  wet_lab_data <- readRDS(paste("mbqc_paper_data_products/",wet_lab_file,sep = ""))
  
  #store data in generic object
  # wet_lab_data <- get(paste("wet_lab_",wet_lab_file,sep = "")) 
  
  #remove non-generic object 
  # rm(list = paste("wet_lab_",wet_lab_file,sep = ""))
  if(sublab_ind){
    suffix <- readRDS(paste(wet_lab,"suffix",sep = "_"))
    suffix <- suffix == as.numeric(wet_sublab) 
   
    
    wet_lab_data <- wet_lab_data[(wet_lab_data$pre.extracted_bool == "FALSE") &
               (wet_lab_data$sequencing_wetlab == wet_lab) &
               !is.na(wet_lab_data$sequencing_wetlab) &
               !is.na(wet_lab_data$pre.extracted_bool),] 
    num_rows <- wet_lab_data %>%
      filter(!is.na(specimen)) %>% 
      filter(sequencing_wetlab == wet_lab) %>% 
      filter(extraction_wetlab == wet_lab | (wet_lab == "HL-F" & extraction_wetlab == "HL-G")) %>%
      dim() %>%
      (function(x) x[1])
    
    if(num_rows != length(suffix)){
      stop("Number of rows in filtered data not equal to length of protocol suffix.")
    } else{
      
      wet_lab_data %<>%
      filter(!is.na(specimen)) %<>% 
      filter(sequencing_wetlab == wet_lab) %<>% 
      filter(extraction_wetlab == wet_lab | (wet_lab == "HL-F" & extraction_wetlab == "HL-G")) %<>%
      filter(suffix)
      }
      
  }
  #filter data 
  wet_lab_data %<>%
    filter(dry_lab %in% included_dry_labs) %<>% #only samples in included dry labs
    filter(pre.extracted_bool == FALSE) %<>% #only non-pre-extracted samples
    filter(sequencing_wetlab == wet_lab) %<>% #only samples sequenced at wet lab in question
    filter((extraction_wetlab == wet_lab)| 
             sequencing_wetlab == "HL-F" & extraction_wetlab == "HL-G") %<>% #only samples extracted at wet lab in question
    filter(!is.na(specimen)) #only valid specimens
    
  
  #get unique specimens
  unique_specimens <- wet_lab_data$specimen %>%
    unique 
  
  #create column for cv fold and for test/train
  wet_lab_data$fold <- NA
  wet_lab_data$test <- NA
  
  #set test/train seed (should be same across all training for same wet lab)
  set.seed(seed)
  
  #split out test data (stratifying by specimen and dry lab to ensure balance)
  ### legacy code archived
  # for(spec in unique_specimens){
  #   for(bl in included_dry_labs){
  #     wet_lab_data$test[
  #       (wet_lab_data$specimen == spec)&
  #         (wet_lab_data$dry_lab == bl)] %<>%
  #       (function(x) sample(rep(c(TRUE,FALSE),ceiling(length(x)/(2))), length(x),replace=FALSE))}
  # }
  ### new code splitting bioinformatics IDS:
  ### first check that no bioinformatics ID is shared across specimens (shouldn't be)
  wet_lab_data %>%
    with(table(Bioinformatics.ID,specimen)) %>%
    apply(1, function(x) sum(x>0)) %>%
    max %>%
    (function(x){ 
      if(x>1){
      stop("Bioinformatics IDs are shared across specimens")} else{
        print("Checked that bioinformatics IDs are not shared across specimens")}
      })
  
  for(spec in unique_specimens){
    unique_bioinformatics.IDs <- wet_lab_data %>%
      filter(specimen == spec) %>%
      with(Bioinformatics.ID) %>%
      unique
    
    #count number of aliquots for specimen
    num_aliquots <- length(unique_bioinformatics.IDs)
    
    #create test or training set assignment for each aliquot
    test_assignments <- sample(rep(c(T,F),
                                   ceiling(num_aliquots/2)),
                               replace = F)[1:num_aliquots]
    
    # add assignments to wet_lab_data
    for(k in 1:num_aliquots){
      wet_lab_data[wet_lab_data$Bioinformatics.ID == unique_bioinformatics.IDs[k],"test"] <- 
        test_assignments[k]
    }
  }
  
  # in the case (HL-C only) where only one sample from a specimen
  # is available, make sure that sample is assigned to the training set
  unique_training <- unique(wet_lab_data[!wet_lab_data$test,"specimen"])
  unique_test <- unique(wet_lab_data[wet_lab_data$test,"specimen"])
  
  if(length(unique_training)< length(unique_test)){
    warning("Not all specimens in test set represented in training set.
Reassigning missing specimens to training set.")
    missing_specimen <- setdiff(unique_test,unique_training)
    wet_lab_data[wet_lab_data$specimen == missing_specimen,"test"] <- FALSE
  }
  
  
  
  print("about to aggregate")
    aggregated_data <- aggregate_level(wet_lab_data,level)
   
    test_data <- aggregated_data[wet_lab_data$test,]
    train_data <- aggregated_data[!wet_lab_data$test,]
    
    if(sublab_ind){
      wet_lab <- paste(wet_lab,wet_sublab,sep="_")
    }

    saveRDS(test_data,paste(wet_lab,"test",level,seed,sep="_"))
    saveRDS(train_data,paste(wet_lab,"train",level,seed,sep="_"))
  
 return("done")

}

### function to transform data to given scale
transform_data <- function(untransformed_data,
                           transformation){
  
  #convert to proportion-scale data if indicated
  if((transformation == "proportion")|
     (transformation == "proportion_diff")|
     transformation == "proportion_diff_scale"|
     transformation == "proportion_ecdf"){
    #convert to relative abundances
    row_totals <- apply(untransformed_data,1,sum)
    for(i in 1:length(row_totals)){
      untransformed_data[i,] <- untransformed_data[i,]/(ifelse(row_totals[i] > 0, row_totals[i], 1))
    }
  }
  
  #convert to presence-absence data if indicated
  if(transformation == "presence"){ 
    untransformed_data %<>% apply(c(1,2), function(x){ 
      if(x>0){return(1)} else{ #return 1 if count > 0
        return(0) #else return 0
      }
    })
  }
  
  #convert to centered-log-ratio transformed data (using pseudocount of 1) if indicated
  if((transformation == "clr_pseudo_1")|
     (transformation == "clr_pseudo_1_diff")|
    ( transformation == "clr_pseudo_1_diff_scale")|
    (transformation == "clr_pseudo_1_ecdf")){ 
    untransformed_data %<>% apply(c(1,2), function(x){
      if(x>0){ return(x)} else{ #if count >0, return count
        return(1) # else return 1 (pseudocount)
      }
    })
    for(i in 1:(nrow(untransformed_data))){
      untransformed_data[i,] %<>% (function(x) log(x) - mean(log(x)))
    }
    if((sum(is.na(untransformed_data))>0) | sum(is.nan(untransformed_data))> 0){
      stop("NANs introduced by transformation")
    }
  }
  
  return("transformed_data" = untransformed_data)
  
}

### function to center data
center_data <- function(transformed_data,
                        specimens){
  
  specimens %<>% as.character()
  
  specs <- unique(specimens)
  
 
  spec_means <- matrix(nrow = length(specs),
                       ncol = ncol(transformed_data))
  
  for(i in 1:length(specs)){
    spec <- specs[i]
    spec_means[i,] <- apply(transformed_data[specimens == spec,,drop = F],2,mean)
  }
  
  spec_mean <- apply(spec_means,2,mean)
  
  for(i in 1:nrow(transformed_data)){
    transformed_data[i,] <- (transformed_data[i,] - spec_mean)
  }
  return(transformed_data)
}

### function to transform to taxon-wise ecdf

ecdf_by_taxon <- function(transformed_data){

  for(j in 1:ncol(transformed_data)){
    transformed_data[,j] %<>% (function(x) sapply(x, function(k) mean(x<=k)))
  }
  return(transformed_data)
}

### function to get between-specimen sds
get_sds <- function(transformed_data,
                    specimens){
  
  specimens %<>% as.character()
  
  specs <- unique(specimens)
  
  
  spec_means <- matrix(nrow = length(specs),
                       ncol = ncol(transformed_data))
  
  for(i in 1:length(specs)){
    spec <- specs[i]
    spec_means[i,] <- apply(transformed_data[specimens == spec,,drop = F],2,mean)
  }
  
  spec_vars <- apply(spec_means,2,function(x) var(x, na.rm= TRUE))
  spec_vars[spec_vars ==0] <- 1 #i.e., leave alone variables that don't vary
  spec_sds <- sqrt(spec_vars)
  return(spec_sds)
  
}

### function to rescale data
rescale_data <- function(transformed_data,
                         spec_sds){
  
  for(i in 1:nrow(transformed_data)){
    transformed_data[i,] <- (transformed_data[i,]/spec_sds)
  }
  return(transformed_data)
}

#function to evaluate xgboost CV performance on training data for one set of parameters
do_one_xgb_cv <- function(wet_lab,
                          nfolds = 10,
                          class_var,
                          eta,
                          max_depth,
                          colsample_bytree,
                          subsample,
                          cv_seed,
                          train_seed,
                          level,
                          transformation = "proportion",
                          dry_lab_subset = NULL){

  
  #print wet lab
  print(wet_lab)
  print(dry_lab_subset)

  
  #load wet lab training data

  train_data <- readRDS(paste(wet_lab,"_train_",level,"_",train_seed,sep = ""))
  
  if(!is.null(dry_lab_subset)){
    train_data %<>% dplyr::filter(dry_lab %in% dry_lab_subset)
  }

  #collapse human specimen types (following Sinha, et al.)
  if(class_var == "specimen_type_collapsed"){
    train_data$specimen_type_collapsed %<>%
      (function(x){ x[x=="Fresh"|x=="Freeze-dried"] <- "Human";
       return(x)})
  }
  
  # return(dim(train_data))
  # print("2")
  #extract classes
  unique_classes <- train_data[,class_var,drop=TRUE] %>% unique
  #make sure classes are in same order
  unique_classes <- unique_classes[order(unique_classes)]
  unique_classes %<>% (function(x) x[!is.na(x)])
  #set cv seed (can vary across parameters)
  set.seed(cv_seed)
  # print("3")
  #construct folds *on aliquots* as xgboost prefers them
  unique_bio_IDs <- train_data$Bioinformatics.ID %>%
    unique
  num_bio_IDs <- length(unique_bio_IDs)
  folds <- sample(rep(1:nfolds,ceiling(num_bio_IDs/nfolds)),replace = F)[1:num_bio_IDs]
  for(k in 1:num_bio_IDs){
    train_data$fold[train_data$Bioinformatics.ID == unique_bio_IDs[k]] <- folds[k]
  }

  cv_folds <- lapply(1:(nfolds), function(x) which(train_data$fold == x))
  
  #format label as xgboost prefers
  train_data <- train_data[order(train_data[,class_var]),]
  train_class_label <- sapply(train_data[,class_var,drop = TRUE], function(x) which(unique_classes == x) -1)
  
  #construct covariate matrix for xgboost
  specimens <- train_data$specimen
  train_data <- train_data %<>%
    select(starts_with("k__")) %<>%
    apply(2,as.numeric)
  
  #transform data to correct scale
  
  train_data %<>% (function(x) transform_data(x,transformation))
  
  if(transformation %in% c("clr_pseudo_1_ecdf",
                           "proportion_ecdf")){
    train_data %<>% ecdf_by_taxon()
  }
  
  
  #center to mean if indicated
  if(transformation == "clr_pseudo_1_diff"|
     transformation == "clr_pseudo_1_diff_scale"|
     transformation == "proportion_diff"|
     transformation == "proportion_diff_scale"){ #subtract specimen-weighted mean
   sds <- get_sds(train_data,specimens = specimens)
   train_data %<>% (function(x) center_data(transformed_data = x, specimens = specimens))

  }
  
  #rescale by weighted sd if indicated
  if(transformation == "clr_pseudo_1_diff_scale"|
     transformation == "proportion_diff_scale"){ #divide by specimen-weighted standard deviations
    train_data %<>% (function(x) rescale_data(transformed_data = x, spec_sds = sds))

  }
  
  
  #convert training data to matrix
  train_data %<>% as.matrix

  # print(dim(train_data))
  # print(length(train_class_label))
  
  #convert to xgb.DMatrix
  train_data %<>%
    xgb.DMatrix(label = train_class_label)
  # print("4")
  all_cv_results <- vector(5,mode = "list")
  counter <- 1
  best_iter <- 0
  best_iters <- numeric(0)
  cv_misclass <- 1
  while((best_iter < 1000) & (counter <= 5) & (cv_misclass > 0)){
    
    print(paste("iteration number",counter))
    #run cross-validation
  all_cv_results[[counter]] <- try(xgb.cv(objective = 'multi:softmax',
                           eval_metric = 'merror',
                           num_class = ifelse(class_var == "specimen",
                                              22,
                                              4),
                           eta = eta,
                           max_depth = max_depth,
                           data = train_data,
                           min_child_weight = 1,
                           lambda = 0,
                           subsample = subsample,
                           gamma = 0,
                           colsample_bytree = colsample_bytree,
                           nrounds = 10000,
                           early_stopping_rounds = 500,
                           folds = cv_folds,
                           nthread = 10,
                           print_every_n = 10)
                           )
  

  #save best iteration
  best_iter <- all_cv_results[[counter]]$best_iteration
  best_iters <- c(best_iters,best_iter)
  
  #collect best cv misclassification
  cv_misclass <- all_cv_results[[counter]]$evaluation_log$test_merror_mean[best_iter]
  
  print(paste("best iterations so far: ", paste(best_iters,collapse = "; "), "\r\n"))
 
  #if best iteration is not changing for at least 4 rounds, cut out of while loop
  
  best_cvs <- all_cv_results[1:counter] %>% sapply(function(x) x$evaluation_log$test_merror_mean[x$best_iteration])
  print(paste("best cv misclassification so far: ", paste(signif(best_cvs,3),collapse = "; "),"\r\n"))

  
  print(paste("best cv misclassification error at step ", best_iter,sep =""))
  
  if(best_iter < 1000){ # if not sufficient number of steps (to meet 1000-step rule of thumb)
    if(cv_misclass >0){ # and if cv misclassification is >0
      if(counter < 5){ # and if we haven't had 10 iterations yet
  
  #set learning rate to smaller value (geometric mean)
  eta <- eta/sqrt(10)
  
  #increment counter

  print(paste("learning rate set to ",eta,sep =""))

      }
    }
  }
  
  counter <- counter + 1
  
  if(cv_misclass == 0){
    print(paste("cv results accepted at # boosting iterations < 1000 as cv misclassification rate (0.00) cannot be improved"))
  }
  }
  
  #save best cv result as main cv result
  # if tie, choose highest number of boosting rounds
  
  for_result <- data_frame(cv = best_cvs,
                           iter = best_iters,
                           trial = 1:(counter-1))
  result_to_select <- for_result %>%
    filter(cv == min(cv)) %>% #only select those trials that had lowest cv misclass error
    filter(iter == max(iter)) %>% #among ties, choose highest number of iterations
    filter(trial == max(trial)) %>% #among any remaining ties, choose smaller eta (later trial)
    select(trial) %>%  #select trial number
    unlist
    
  cv_results <- all_cv_results[[result_to_select]] #select best trial 

  return(list("wet_lab" = wet_lab,
              "level" = level,
              "subsample" = subsample,
              "colsample_bytree" = colsample_bytree,
              "eta" = cv_results$params$eta,
              "max_depth" = max_depth,
              "cv_results" = cv_results,
              "all_cv_results" = all_cv_results,
              "dry_lab_subset" = dry_lab_subset))
  }

#function to evaluate glmnet performance on training data for one set of parameters
do_one_glmnet_cv <- function(wet_lab,
                          nfolds = 10,
                          class_var,
                          alpha,
                          cv_seed,
                          train_seed,
                          level, 
                          transformation = "proportion",
                          included_dry_labs,
                          dry_lab_subset = NULL){

  
  #print wet lab
  print(wet_lab)
  print(dry_lab_subset)
  
  
  #load wet lab training data
  
  train_data <- readRDS(paste(wet_lab,"_train_",level,"_",train_seed,sep = ""))
  
  if(!is.null(dry_lab_subset)){
    train_data %<>% dplyr::filter(dry_lab %in% dry_lab_subset)
  }
  
  #collapse human specimen types (following Sinha, et al.)
  if(class_var == "specimen_type_collapsed"){
    train_data$specimen_type_collapsed %<>%
      (function(x){ x[x=="Fresh"|x=="Freeze-dried"] <- "Human";
      return(x)})
  }
  

  
  # return(dim(train_data))
  # print("2")
  #extract classes
  unique_classes <- train_data[,class_var,drop=TRUE] %>% unique
  # print(unique_classes)
  unique_classes %<>% (function(x) x[!is.na(x)])
  #make sure classes are in same order
  unique_classes <- unique_classes[order(unique_classes)]
  #set cv seed (can vary across parameters)
  set.seed(cv_seed)
  # print("3")
  #construct folds *on aliquots* as xgboost prefers them
  unique_bio_IDs <- train_data$Bioinformatics.ID %>%
    unique
  if(nfolds != "loocv"){
    
  
  num_bio_IDs <- length(unique_bio_IDs)
  folds <- sample(rep(1:nfolds,ceiling(num_bio_IDs/nfolds)),replace = F)[1:num_bio_IDs]
  } else{ # if loocv, each bioinformatics ID gets its own fold
    num_bio_IDs <- length(unique_bio_IDs)
  folds <- 1:num_bio_IDs
  nfolds <- length(folds)
    
  }
  
  for(k in 1:num_bio_IDs){
    train_data$fold[train_data$Bioinformatics.ID == unique_bio_IDs[k]] <- folds[k]
  }
  
  cv_folds <- train_data$fold
  
  train_data[,class_var] %<>% as.character 
  #format label in analogy with what xgboost prefers (but as factor for glmnet)
  train_data <- train_data[order(train_data[,class_var]),]
  train_class_label <- sapply(train_data[,class_var,drop = TRUE], function(x) which(unique_classes == x) -1) 
  
  # #print(train_data[,100] %>%summary)
  
  #construct covariate matrix for xgboost
  specimens <- train_data$specimen
  train_data <- train_data %<>%
    select(starts_with("k__")) %<>%
    apply(2,as.numeric)
  
  #transform data to correct scale
  train_data %<>% (function(x) transform_data(x,transformation))
  
  if(transformation %in% c("clr_pseudo_1_ecdf",
                           "proportion_ecdf")){
    train_data %<>% ecdf_by_taxon()
  }
  
  #center to mean if indicated
  if(transformation == "clr_pseudo_1_diff"|
     transformation == "clr_pseudo_1_diff_scale"|
     transformation == "proportion_diff"|
     transformation == "proportion_diff_scale"){ #subtract specimen-weighted mean
    sds <- get_sds(train_data,specimens = specimens)
    train_data %<>% (function(x) center_data(x, specimens))
  }
  
  #rescale by weighted sd if indicated
  if(transformation == "clr_pseudo_1_diff_scale"|
     transformation == "proportion_diff_scale"){ #divide by specimen-weighted standard deviations
    train_data %<>% (function(x) rescale_data(x, spec_sds = sds))
  }

  #convert to matrix
  train_data %<>% as.matrix
  
  # print(dim(train_data))
  # print(length(train_class_label))


  # 
 
  
  # 
  # 
  unique_classes <- train_class_label %>% unique %>% (function(x) x[order(x)])
  # get max number of 
  max_fold_counts <- unique_classes %>%
    sapply(function(x) sapply(1:nfolds, function(y) sum(train_class_label[cv_folds ==y]==x))) %>%
    apply(2,max)
  class_counts <- unique_classes %>%
    sapply(function(x) sum(train_class_label == x))
  
  min_train_counts <- class_counts - max_fold_counts
  # do any classes have too few observations for glmnet to work?
  class_insuff_data <- unique_classes[min_train_counts < 2]
  # 
  #if so eliminate these classes
  if(length(class_insuff_data)>0){
  train_data <- train_data[!(train_class_label %in% class_insuff_data),]
  cv_folds <- cv_folds[!(train_class_label %in% class_insuff_data)]


  train_class_label <- train_class_label[!(train_class_label %in% class_insuff_data)]

  train_class_label %<>% as.factor
  levels(train_class_label) <- 0:(length(unique_classes) - length(class_insuff_data) -1)
  train_class_label %<>% as.numeric
  }
  # lambdas <-  exp(seq(-10,10,by=.1))
  
  cv_results <- try(cv.glmnet(y = matrix(as.factor(train_class_label), ncol = 1),
                          x = train_data,
                          type.measure = "class",
                          foldid = cv_folds,
                          family = 'multinomial',
                          alpha = alpha))
  
   #assign NA if no convergence
  if(inherits(cv_results,"try-error")){
    cv_results <- NA
  }
  
  # plot(cv_results)

  
  
  return(list("wet_lab" = wet_lab,
              "alpha" = alpha,
              "cv_results" = cv_results,
              "excluded_classes" = class_insuff_data,
              "dry_lab_subset" = dry_lab_subset))
}

### calculate sample-centered centered-log-ratio-transformed data
clr_diff <- function(count_data){
  
  test_data <- count_data
  test_data[,colnames(test_data) %>% sapply(function(x) substr(x,1,3)) %>% (function(x) x== "k__")] %<>%
    apply(2, as.numeric) %<>%
    apply(c(1,2), function(x) if(x == 0){ return(1)} else{return(x)}) %<>% #add pseudocount
    apply(c(1,2), log) %<>% #log transform
    apply(1, function(x) x- mean(x)) #center
  
  
  ncols <- ncol(test_data %>% select(starts_with('k__')))
  specs <- test_data$specimen %>% unique
  spec_means <- matrix(nrow = length(specs),ncol = ncols)
  
  #get mean for each specimen
  for(q in 1:length(specs)){
    spec <- specs[q]
  spec_means[q,] <- (test_data %>%
      filter(specimen == spec) %>%
      select(starts_with("k__")) %>%
      apply(2,mean))
        
  }
  
  spec_mean <- apply(spec_means,2,mean)
  
  for(i in 1:nrow(test_data)){
  test_data[i,colnames(test_data) %>% sapply(function(x) substr(x,1,3)) %>% (function(x) x== "k__")] <- 
    (test_data[i,colnames(test_data) %>% sapply(function(x) substr(x,1,3)) %>% (function(x) x== "k__")] - 
       spec_mean
   )
  }
  
  return(test_data)
}

#function to train glmnet on training set using selected parameters 
#and predict on test sets
do_one_glmnet_train <- function(wet_lab,
                             class_var,
                             alpha,
                             lambda,
                             train_seed,
                             level,
                             included_wet_labs,
                             included_dry_labs,
                             transformation = "proportion",
                             dry_lab_subset = NULL){
  
  print(paste("wet lab: ", wet_lab))
  print(paste("class_var: ", class_var))
  print(paste("alpha: ", alpha))
  print(paste("train_seed: ", train_seed))
  print(paste("taxon: ", level))

  
  #print wet lab
  print(wet_lab)
  if(!is.null(dry_lab_subset)){
  print(paste(c("dry_lab_subset",dry_lab_subset),
              sep = " ", collapse = " "))
    }
  
  
  #load wet lab training data
  
  train_data <- readRDS(paste(wet_lab,"_train_",level,"_",train_seed,sep = ""))
  
  if(!is.null(dry_lab_subset)){
    train_data %<>% dplyr::filter(dry_lab %in% dry_lab_subset)
  }
  
  #collapse human specimen types (following Sinha, et al.)
  if(class_var == "specimen_type_collapsed"){
    train_data$specimen_type_collapsed %<>%
      (function(x){ x[x=="Fresh"|x=="Freeze-dried"] <- "Human";
      return(x)})
  }
  
  # return(dim(train_data))
  # print("2")
  #extract classes
  unique_classes <- train_data[,class_var,drop=TRUE] %>% unique
  # print(unique_classes)
  unique_classes %<>% (function(x) x[!is.na(x)]) %>% unlist 
  #make sure classes are in same order
  unique_classes <- unique_classes[order(unique_classes)]
  #save training classes for use with test data
  unique_classes_train <- unique_classes


  print("about to characterize")
  train_data %<>% as.data.frame()

  for_order <- train_data[,class_var] %>% unlist %>% as.vector %>% as.character

  print("characterized")
  #format label in analogy with what xgboost prefers (but as factor for glmnet)
  train_data <- train_data[order(for_order),]
  train_class_label <- sapply(train_data[,class_var,drop = TRUE], function(x) which(unique_classes == x) -1) 
  
  #prep data for glmnet
  specimens <- train_data$specimen
  train_data <- train_data %<>%
    select(starts_with("k__")) %<>%
    apply(2,as.numeric)
  #transform data to correct scale
  
  train_data %<>% (function(x) transform_data(x,transformation))
  
  if(transformation %in% c("clr_pseudo_1_ecdf",
                           "proportion_ecdf")){
    train_data %<>% ecdf_by_taxon()
  }
  
  #center to mean if indicated
  if(transformation == "clr_pseudo_1_diff"|
     transformation == "clr_pseudo_1_diff_scale"|
     transformation == "proportion_diff"|
     transformation == "proportion_diff_scale"){
    #subtract specimen-weighted mean
    sds <- get_sds(train_data,specimens)
    train_data %<>% (function(x) center_data(x, specimens))
  }
  
  #rescale by weighted sd if indicated
  if(transformation == "clr_pseudo_1_diff_scale"|
     transformation == "proportion_diff_scale"){ #divide by specimen-weighted standard deviations
     train_data %<>% (function(x) rescale_data(x, spec_sds = sds))
  }
  
  
  #convert to matrix
  train_data %<>% as.matrix
  
  unique_classes <- train_class_label %>% unique %>% (function(x) x[order(x)])
  # get max number of 

  class_counts <- unique_classes %>%
    sapply(function(x) sum(train_class_label == x))
  
  min_train_counts <- class_counts
  # do any classes have too few observations for glmnet to work?
  class_insuff_data <- unique_classes[min_train_counts < 2]
  # 
  #if so eliminate these classes
  if(length(class_insuff_data)>0){
    train_data <- train_data[!(train_class_label %in% class_insuff_data),]
    
    
    train_class_label <- train_class_label[!(train_class_label %in% class_insuff_data)]
    
    train_class_label %<>% as.factor
    levels(train_class_label) <- 0:(length(unique_classes) - length(class_insuff_data) -1)
    train_class_label %<>% as.numeric
  }

  print(paste("lambda =", lambda))
  print(paste("alpha = ", alpha))
  
  ##print(train_data[,100] %>%summary)
  
  glmnet_results <- try(glmnet(y = matrix(as.factor(train_class_label), ncol = 1),
                              x = train_data,
                              family = 'multinomial',
                              alpha = alpha,
                              lambda = exp(seq(9,-1,by=-.1))*lambda))
  


  
  glmnet_predictions <- data.frame(
    wet_lab = character(0),
    true_labels = character(0),
    pred_labels = character(0)
  )
  
  glmnet_misclass <- data.frame(
    wet_lab = character(0),
    misclass = numeric(0)
  )
  
  for(wl in included_wet_labs){
    print(paste("predicting on", wl,sep = " "))
    test_data <- readRDS(paste(wl,"_test_",level,"_",train_seed,sep = ""))
    
    
    if(!is.null(dry_lab_subset)){
      test_data %<>% dplyr::filter(dry_lab %in% dry_lab_subset)
    }
    
    if(class_var == "specimen_type_collapsed"){
      test_data$specimen_type_collapsed %<>%
        (function(x){ x[x=="Fresh"|x=="Freeze-dried"] <- "Human";
        return(x)})
    }
    
    # return(dim(train_data))
    # print("2")
    #extract classes
    unique_classes <- test_data[,class_var,drop=TRUE] %>% unique
    # print(unique_classes)
    unique_classes %<>% (function(x) x[!is.na(x)]) %>% unlist 
    #make sure classes are in same order
    unique_classes <- unique_classes[order(unique_classes)]
    
    
  
    test_data %<>% as.data.frame()
    
    for_order <- test_data[,class_var] %>% unlist %>% as.vector %>% as.character
    
    #format label in analogy with what xgboost prefers (but as factor for glmnet)
    test_data <- test_data[order(for_order),]
    test_class_label <- sapply(test_data[,class_var,drop = TRUE], function(x) which(unique_classes_train == x) -1) 

    #construct covariate matrix for xgboost
    specimens_test <- test_data$specimen
    test_data <- test_data %<>%
      select(starts_with("k__")) %<>%
      apply(2,as.numeric)
    
    #transform data to correct scale
    test_data %<>% (function(x) transform_data(x,transformation))
    
    if(transformation %in% c("clr_pseudo_1_ecdf",
                             "proportion_ecdf")){
      test_data %<>% ecdf_by_taxon()
    }
    
    #center to mean if indicated
    if(transformation == "clr_pseudo_1_diff"|
       transformation == "clr_pseudo_1_diff_scale"|
       transformation == "proportion_diff"|
       transformation == "proportion_diff_scale"){ #subtract specimen-weighted mean
      test_sds <- get_sds(test_data, specimens_test)
      test_data %<>% (function(x) center_data(x, specimens_test))
    }
    
    #rescale by weighted sd if indicated
    if(transformation == "clr_pseudo_1_diff_scale"|
       transformation == "proportion_diff_scale"){ #divide by specimen-weighted standard deviations
       test_data %<>% (function(x) rescale_data(x, test_sds))
    }
    
    
    #convert to matrix 
    test_data %<>% as.matrix
    
    # print(dim(train_data))
    # print(length(train_class_label))
    preds <- predict(glmnet_results,newx = test_data,
                     s = lambda,
                     type = "class") %>% as.vector()
    
    
    glmnet_predictions_new <- data.frame(wet_lab = rep(wl,length(preds)),
                                         true_labels = test_class_label,
                                         pred_labels = preds)
    
    glmnet_predictions <- rbind(glmnet_predictions,glmnet_predictions_new)
    # 
    
    glmnet_misclass_new <- data.frame(wet_lab = wl,
                                      misclass = mean(preds != test_class_label))
    
    glmnet_misclass <- rbind(glmnet_misclass,glmnet_misclass_new)
    
    # 
    
  }

  return(list("wet_lab" = wet_lab,
              "class_var" = class_var,
              "taxon" = level,
              "glmnet_results" = glmnet_results,
              "glmnet_misclass" = glmnet_misclass,
              "glmnet_preds" = glmnet_predictions,
              "excluded_classes" = class_insuff_data,
              "dry_lab_subset" = dry_lab_subset))
}

#function to collect best performing glmnet parameters,
# train on training sets with them, and predict on all test sets
train_glmnet <- function(glmnet_params, 
                         included_wet_labs, 
                         transformation = "proportion",
                         queue = "students.q",
                         dry_lab_subset = NULL){
  
  if(is.null(dry_lab_subset)){
    cv_dir <- paste("registry_paper_glmnet",transformation, sep = "_")
    train_dir <- paste("registry_glmnet_train",transformation,sep = "_")
  } else{
    cv_dir <- paste("registry_paper_glmnet",transformation,
                    dry_lab_subset, sep = "_")
    train_dir <- paste("registry_glmnet_train",transformation,
                       dry_lab_subset, sep = "_")
  }

  
  if((dir.exists(cv_dir)) &
     (!dir.exists(train_dir))){
    
  if(transformation == "lr_pseudo_1"|transformation == "lr_pseudo_halfsmallest"){
    # glmnet_params <- expand.grid(included_wet_labs,
    #                              alpha = seq(0,1, by = .05),
    #                              c("specimen","specimen_type_collapsed"),
    #                              c(
    #                                "phylum",
    #                                "order",
    #                                "class",
    #                                "family"))
  }
  nglmnet_results <- dim(glmnet_params)[1]
  glmnet_results <- data.frame(wet_lab = character(nglmnet_results ),
                               alpha = numeric(nglmnet_results),
                               lambda = numeric(nglmnet_results),
                               cv = numeric(nglmnet_results ),
                               cv_sd = numeric(nglmnet_results),
                               class_var = glmnet_params[1:nglmnet_results,3],
                               taxon = as.character(glmnet_params[1:nglmnet_results,4]))
  glmnet_results$wet_lab %<>% as.character()
  glmnet_results$alpha %<>% as.numeric()
  glmnet_results$lambda %<>% as.numeric()
  glmnet_results$class_var %<>% as.character()
  glmnet_results$taxon %<>% as.character

  for(i in 1:nglmnet_results ){
    print(paste("glmnet result ", i, " of ", nglmnet_results))
    job <- try(readRDS(paste(cv_dir,"/results/",i,".rds",sep = "")))
    if(is.list(job)){
      if(is.list(job$cv_results)){
        glmnet_results$wet_lab[i] <- as.character(job$wet_lab)
        glmnet_results$alpha[i] <- job$alpha
        glmnet_results$lambda[i] <- job$cv_results$lambda.min
        glmnet_results$cv[i] <-  (job$cv_results$cvm %>% min)
        glmnet_results$cv_sd[i] <- job$cv_results$cvsd[which.min(job$cv_results$cvm)]

      }} else{
        glmnet_results[i,c("lambda","cv","cv_sd")] <- NA
        glmnet_results$wet_lab[i] <- as.character(job$wet_lab)
        glmnet_results$alpha[i] <- job$alpha
      }
    if(glmnet_results[i,"wet_lab"] == ""){
      glmnet_results[i,"wet_lab"] <- as.character(glmnet_params[i,1])
      glmnet_results[i,"cv"] <- 100
    }
    print(glmnet_results[i,])
  }
  
  
  glmnet_results$wet_lab %<>% as.character %<>% as.factor
  
  if(is.null(dry_lab_subset)){
  saveRDS(glmnet_results, file = paste("glmnet_cv_results", transformation, sep = "_",collapse = "_"))
  } else{
    saveRDS(glmnet_results, file = paste("glmnet_cv_results", transformation,
                                         dry_lab_subset,
                                         sep = "_",collapse = "_"))
  }
  
  
  ### summarize 10-fold cv results
  best_glmnet <-
    glmnet_results %>%
    group_by(class_var) %>%
    group_by(wet_lab, .add = TRUE) %>%
    group_by(taxon,.add = TRUE) %>%
    filter(cv == min(cv)) %>%
    summarize(cv = min(cv),
              alpha = paste(unique(alpha[cv == min(cv)]),collapse = "; "))
  
 

  
  
  alpha_list <- best_glmnet$alpha %>%
    sapply(function(x) #split out multiple values of alpha into character vector
      base::strsplit(x, split = "; ", fixed = TRUE)) %>%
    lapply(as.numeric)  #make numeric again
  
  names(alpha_list) <- NULL
   
  best_glmnet$alpha <- alpha_list %>%
    sapply(function(x) x[abs(x-.49) == min(abs(x-.49))]) %>%  #among ties, pull out
    unlist
  #alpha closest to .5,
  #resolving ties downward
  #find job ID
  glmnet_results$job_id <- 1:nglmnet_results
  best_glmnet$job_id <-
    best_glmnet %>%
    apply(1, function(x){
      ((glmnet_results$wet_lab == unlist(x["wet_lab"])) &
         (glmnet_results$taxon == unlist(x["taxon"])) &
         (glmnet_results$class_var == unlist(x["class_var"])) &
         (abs(glmnet_results$alpha - as.numeric(unlist(x["alpha"])))<(1e-4)))   %>%
        which
    }) %>%
    unlist
  
  #peel out results from glmnet_results now and save in best_glmnet
  best_glmnet <-
    glmnet_results[best_glmnet$job_id,]
  
  print(best_glmnet)
  
  best_glmnet %<>%
    filter(wet_lab != "")
  
  print(best_glmnet)
  
   #run only if results not already calculated
  if(is.null(dry_lab_subset)){
    
    glmnet_train_reg = makeRegistry(file.dir = paste("registry_glmnet_train",transformation,sep = "_"),
                                    packages = c("magrittr",
                                                 "xgboost",
                                                 "tidyr",
                                                 "dplyr",
                                                 "RDS",
                                                 "data.table",
                                                 "glmnet"),
                                    source = "all_functions.R")
  } else{
    glmnet_train_reg = makeRegistry(file.dir = paste("registry_glmnet_train",
                                                     transformation,
                                                     dry_lab_subset, sep = "_"),
                                    packages = c("magrittr",
                                                 "xgboost",
                                                 "tidyr",
                                                 "dplyr",
                                                 "RDS",
                                                 "data.table",
                                                 "glmnet"),
                                    source = "all_functions.R")
    }
    
    options(batchtools.progress = FALSE)
    print(best_glmnet)
    print(dim(best_glmnet))
    
    batchMap(fun = do_one_glmnet_train,
             wet_lab = best_glmnet$wet_lab,
             alpha = best_glmnet$alpha,
             lambda = best_glmnet$lambda,
             class_var = best_glmnet$class_var,
             level = best_glmnet$taxon,
             train_seed = train_seed,
             reg = glmnet_train_reg,
             more.args = list(included_wet_labs = included_wet_labs,
                              transformation = transformation,
                              included_dry_labs = included_dry_labs,
                              dry_lab_subset = dry_lab_subset))
    
    submitJobs(reg = glmnet_train_reg,
               resources = list(memory = "25G",
                                queue = queue))
    
    
    waitForJobs()
    
    #delete (large) directory containing cross-validation results
    print(paste("Deleting ", "registry_paper_glmnet_",transformation,
                " to free memory",sep = ""))
    
    if(is.null(dry_lab_subset)){

    unlink(paste("registry_paper_glmnet_",transformation,sep = ""), recursive = TRUE)
    } else{
        
      unlink(paste("registry_paper_glmnet_",transformation,
                   dry_lab_subset,
                   sep = ""), recursive = TRUE)
      }
    
  }
  
}


do_one_xgb_train <- function(wet_lab,
                          class_var,
                          eta,
                          max_depth,
                          colsample_bytree,
                          subsample,
                          train_seed,
                          cv_seed,
                          level,
                          best_iteration,
                          included_wet_labs,
                          dry_lab_subset,
                          transformation = "proportion"){
  
  #load libraries
  library(xgboost)
  library(magrittr)
  library(tidyr)
  library(dplyr)
  library(RDS)
  library(data.table)
  
  #print wet lab
  print(wet_lab)
  if(!is.null(dry_lab_subset)){
    print(paste(c("dry_lab_subset",dry_lab_subset),
                sep = " ", collapse = " "))
  }
  
  
  #load wet lab training data
  
  train_data <- readRDS(paste(wet_lab,"_train_",level,"_",train_seed,sep = ""))
  
  if(!is.null(dry_lab_subset)){
    train_data %<>% dplyr::filter(dry_lab %in% dry_lab_subset)
  }

  class_var %<>% as.character
  print(class_var)
  print(train_data[1:4,class_var])
  #collapse human specimen types (following Sinha, et al.)
  if(class_var == "specimen_type_collapsed"){
    train_data$specimen_type_collapsed %<>%
      (function(x){ x[x=="Fresh"|x=="Freeze-dried"] <- "Human";
      return(x)})
  }
  
  # return(dim(train_data))
  # print("2")
  #extract classes
  unique_classes <- train_data[,class_var,drop=TRUE] %>% unique
  #make sure classes are in same order
  unique_classes <- unique_classes[order(unique_classes)]
  # print(unique_classes)
  unique_classes %<>% (function(x) x[!is.na(x)])
  #set cv seed (can vary across parameters)
  set.seed(cv_seed)
  # print("3")
  #construct folds as xgboost prefers them
  # train_data$fold %<>%
    # (function(x) sample(rep(1:nfolds, ceiling(length(x)/nfolds)),length(x),replace = FALSE))
  
  # cv_folds <- lapply(1:(nfolds), function(x) which(train_data$fold == x))
  
  #format label as xgboost prefers
  train_data <- train_data[order(train_data[,class_var]),]
  train_class_label <- sapply(train_data[,class_var,drop = TRUE], function(x) which(unique_classes == x) -1)
  
  #construct covariate matrix for xgboost
  specimens <- train_data$specimen
  train_data <- train_data %<>%
    select(starts_with("k__")) %<>%
    apply(2,as.numeric)
  
  
  #convert to proportion-scale data if indicated
  train_data %<>% (function(x) transform_data(x,transformation))
  
  if(transformation %in% c("clr_pseudo_1_ecdf",
                           "proportion_ecdf")){
    train_data %<>% ecdf_by_taxon()
  }
  
  #center to mean if indicated
  if(transformation == "clr_pseudo_1_diff"|
     transformation == "clr_pseudo_1_diff_scale"|
     transformation == "proportion_diff"|
     transformation == "proportion_diff_scale"){ #subtract specimen-weighted mean
    sds <- get_sds(train_data,specimens)
    train_data %<>% (function(x) center_data(x, specimens))
  }
  
  #rescale by weighted sd if indicated
  if(transformation == "clr_pseudo_1_diff_scale"|
     transformation == "proportion_diff_scale"){ #divide by specimen-weighted standard deviations
    train_data %<>% (function(x) rescale_data(x, sds))
  }
  
  
  #convert training data to matrix
  train_data %<>% as.matrix
  
  # print(dim(train_data))
  # print(length(train_class_label))
  
  #convert to xgb.DMatrix
  train_data %<>%
    xgb.DMatrix(label = train_class_label)
  
  print(summary(train_class_label))
  # print("4")
  
  xgb_results <- try(xgb.train(objective = 'multi:softmax',
                           num_class = ifelse(class_var == "specimen",
                                              22,
                                              4),
                           eta = eta,
                           max_depth = max_depth,
                           data = train_data,
                           min_child_weight = 1,
                           lambda = 0,
                           subsample = subsample,
                           gamma = 0,
                           colsample_bytree = colsample_bytree,
                           nrounds = best_iteration,
                           nthread = 10))
  

  
  
  
  xgb_predictions <- data.frame(
    wet_lab = character(0),
    true_labels = character(0),
    pred_labels = character(0)
  )
  
  xgb_misclass <- data.frame(
    wet_lab = character(0),
    misclass = numeric(0)
  )
  
  for(wl in included_wet_labs){
    print(paste("predicting on", wl,sep = " "))
    test_data <- readRDS(paste(wl,"_test_",level,"_",train_seed,sep = ""))
    
    if(!is.null(dry_lab_subset)){
      test_data %<>% dplyr::filter(dry_lab %in% dry_lab_subset)
    }
    
    if(class_var == "specimen_type_collapsed"){
      test_data$specimen_type_collapsed %<>%
        (function(x){ x[x=="Fresh"|x=="Freeze-dried"] <- "Human";
        return(x)})
    }
    
    # return(dim(train_data))
    # print("2")
    #extract classes
    # unique_classes <- test_data[,class_var,drop=TRUE] %>% unique
    # print(unique_classes)
    # unique_classes %<>% (function(x) x[!is.na(x)]) %>% unlist 
    #make sure classes are in same order
    # unique_classes <- unique_classes[order(unique_classes)]
    #set cv seed (can vary across parameters)
    
    
    print("about to characterize")
    test_data %<>% as.data.frame()
    
    for_order <- test_data[,class_var] %>% unlist %>% as.vector %>% as.character
    
    print("characterized")
    test_data <- test_data[order(test_data[,class_var]),]
    test_class_label <- sapply(test_data[,class_var,drop = TRUE], function(x) which(unique_classes == x) -1)
    
    # if(transformation == "clr_pseudo_1_diff"){
    #   test_data %<>% clr_diff
    #   test_data %>%
    #     select(starts_with("k__")) %>%
    #     apply(2, mean) %>%
    #     (function(x) max(abs(x))) %>%
    #     (function(x) paste("max abs mean", signif(x,3)))
    # }
    
    
    #construct covariate matrix for xgboost
    specimens_test <- test_data$specimen
    dry_labs_test <- test_data$dry_lab
    test_data <- test_data %<>%
      select(starts_with("k__")) %<>%
      apply(2,as.numeric)
    
    
    #transform data to correct scale
    test_data %<>% (function(x) transform_data(x,transformation))
    
    if(transformation %in% c("clr_pseudo_1_ecdf",
                             "proportion_ecdf")){
      test_data %<>% ecdf_by_taxon()
    }
    
    
    
    #center to mean if indicated
    if(transformation == "clr_pseudo_1_diff"|
       transformation == "clr_pseudo_1_diff_scale"|
       transformation == "proportion_diff"|
       transformation == "proportion_diff_scale"){ #subtract specimen-weighted mean
      test_sds <- get_sds(test_data,specimens_test)
      test_data %<>% (function(x) center_data(x, specimens_test))
    }
    
    #rescale by weighted sd if indicated
    if(transformation == "clr_pseudo_1_diff_scale"|
       transformation == "proportion_diff_scale"){ #divide by specimen-weighted standard deviations
      test_data %<>% (function(x) rescale_data(x, test_sds))
    }
    
    
    
    #convert to matrix
    test_data %<>% as.matrix
    
    # print(dim(train_data))
    # print(length(train_class_label))
    preds <- predict(xgb_results,newdata = test_data)
    
    
    xgb_predictions_new <- data.frame(wet_lab = rep(wl,length(preds)),
                                         true_labels = test_class_label,
                                         pred_labels = preds)
    
   xgb_predictions <- rbind(xgb_predictions,xgb_predictions_new)
    # 
    
    xgb_misclass_new <- data.frame(wet_lab = wl,
                                      misclass = mean(preds != test_class_label))
    
    xgb_misclass <- rbind(xgb_misclass,xgb_misclass_new)
    
    # 
    
  }
  
  
  return(list("wet_lab" = wet_lab,
              "class_var" = class_var,
              "taxon" = level,
              "xgb_results" = xgb_results,
              "xgb_misclass" = xgb_misclass,
              "xgb_preds" = xgb_predictions,
              "dry_lab_subset" = dry_lab_subset))
}

  


#function to collect best performing xgboost parameters,
# train on training sets with them, and predict on all test sets
train_xgb <- function(xgb_params,
                      included_wet_labs,
                      transformation = "proportion",
                      cv_seed = 5466,
                      dry_lab_subset,
                      resources){
  # check that training registry does not already exist
  # and that ??? registry does exist
  if(is.null(dry_lab_subset)){
    cv_dir <- paste("registry_paper_xgb",transformation, sep = "_")
    train_dir <-paste("registry_xgb_train",transformation,sep = "_")
  } else{
    cv_dir <- paste("registry_paper_xgb",transformation,
                    dry_lab_subset,sep = "_")
    
    train_dir <- paste("registry_xgb_train",transformation,
                       dry_lab_subset, sep = "_")
  }
  
  print(dry_lab_subset)
  
  
  if((dir.exists(cv_dir)) &
     (!dir.exists(train_dir))){
    
  nxgb_results <- dim(xgb_params)[1]
  xgb_results <- data.frame(wet_lab = character(nxgb_results ),
                            eta = numeric(nxgb_results),
                            max_depth = numeric(nxgb_results ),
                            colsample_bytree =  xgb_params[1:nxgb_results,5],
                            cv = numeric(nxgb_results ),
                            iter = numeric(nxgb_results ),
                            subsample = xgb_params[1:nxgb_results,4],
                            class_var = xgb_params[1:nxgb_results,6],
                            taxon = factor(nxgb_results,levels = unique(xgb_params[,7])),
                            cv_sd <- numeric(nxgb_results))
  xgb_results$wet_lab %<>%
    as.character()
  factor(levels =paste("HL-",LETTERS[1:14],sep="") )
  xgb_results$eta %<>% as.numeric()
  xgb_results$colsample_bytree %<>% as.numeric()
  xgb_results$cv %<>% as.numeric()
  xgb_results$subsample %<>% as.numeric()
  for(i in 1:nxgb_results ){
    job <- try(readRDS(paste(cv_dir,"/results/",i,".rds",sep = "")))
    if(is.list(job)){
      xgb_results$wet_lab[i] <- job$wet_lab
      xgb_results$taxon[i] <- job$level
      xgb_results$eta[i] <- job$cv_results$params$eta
      xgb_results$max_depth[i] <-  job$cv_results$params$max_depth
      xgb_results$cv[i] <- job$cv_results$evaluation_log$test_merror_mean[
        job$cv_results$best_iteration]
      xgb_results$iter[i] <- (job$cv_results$evaluation_log$iter[
        job$cv_results$evaluation_log$test_merror_mean == xgb_results$cv[i]] %>%
          median %>%
          floor
      )
      xgb_results$colsample_bytree[i] <- job$cv_results$params$colsample_bytree
      xgb_results$cv_sd[i] <- job$cv_results$evaluation_log$test_merror_std[
        job$cv_results$best_iteration]
    } else{
      xgb_results[i,] <- NA
    }
  }

  
  # get best 10-fold CV results on training data by lab
  best_xgb <-  xgb_results %>%
    # (function(x) x[complete.cases(x),]) %>%
    group_by(wet_lab) %>% #for each combination of wet lab
    group_by(taxon, add = TRUE) %>% #taxon
    group_by(class_var, add = TRUE) %>% #and classification task
    summarize(cv = min(cv),
              colsample_bytree = colsample_bytree[which.min(cv)],
              subsample = subsample[which.min(cv)],
              iter = iter[which.min(cv)],
              eta = eta[which.min(cv)])
  

  # resolve conflicts if there are ties for best cv misclassification error
  for(k in 1:nrow(best_xgb)){
    best_xgb[k, c("iter", "eta","subsample","colsample_bytree")] <- 
      xgb_results %>%
      filter(wet_lab == best_xgb$wet_lab[k]) %>% #look at results only for current wet lab
      filter(taxon == best_xgb$taxon[k]) %>% #taxon
      filter(class_var == best_xgb$class_var[k]) %>% #and classification task
      filter(cv == best_xgb$cv[k]) %>% #restrict to those that made the bar
      #order by largest iteration, then by smallest eta, subsample, colsample by tree
      (function(x) x[order(x$iter, -x$eta, -x$subsample, -x$colsample_bytree,decreasing = TRUE),,drop=FALSE]) %>%
      #pull out parameters
      (function(x) x[1,c("iter", "eta","subsample","colsample_bytree")])
      
  }
  
  
  
  #find job ID
  xgb_results$job_id <- 1:nxgb_results
  best_xgb$job_id <-
    best_xgb %>%
    apply(1, function(x){
      ((xgb_results$wet_lab == unlist(x["wet_lab"])) &
         (xgb_results$taxon == unlist(x["taxon"])) &
         (xgb_results$class_var == unlist(x["class_var"])) &
         (abs(xgb_results$iter - as.numeric(unlist(x["iter"])))<(1e-4)) &
         (abs(xgb_results$colsample_bytree - as.numeric(unlist(x["colsample_bytree"])))<(1e-4)) &
         (abs(xgb_results$subsample - as.numeric(unlist(x["subsample"])))<(1e-4)))   %>%
        which
    }) %>%
    unlist
  
  #just peel out results fromm xgb_results now and save in best_xgb
  best_xgb <-
    xgb_results[best_xgb$job_id,]
  

  
  
  if(!dir.exists(train_dir)){ #run only if results not already calculated
    
    xgb_train_reg = makeRegistry(file.dir = train_dir,
                                 conf.file = "batchtools_multi.conf.R",
                                 packages = c("magrittr",
                                              "xgboost",
                                              "tidyr",
                                              "dplyr",
                                              "RDS",
                                              "data.table",
                                              "xgboost"),
                                 source = "all_functions.R"
    )
    
    options(batchtools.progress = TRUE)
    
    
    batchMap(fun = do_one_xgb_train,
             wet_lab = best_xgb$wet_lab,
             class_var = best_xgb$class_var,
             eta = best_xgb$eta,
             colsample_bytree = best_xgb$colsample_bytree,
             subsample = best_xgb$subsample,
             level = best_xgb$taxon,
             best_iteration = best_xgb$iter,
             reg = xgb_train_reg,
             more.args = list(included_wet_labs = included_wet_labs,
                              transformation = transformation,
                              train_seed = train_seed,
                              max_depth = 6,
                              cv_seed = cv_seed,
                              dry_lab_subset = dry_lab_subset))
    
    
                         
    submitJobs(reg = xgb_train_reg,
               resources = resources)
    
    
    waitForJobs()
    
  }
  }
}




# function to evaluate cross-validated performance of
# glmnet parameters in parallel via batchtools

glmnet_cross_validate <-  function(included_wet_labs,
                                   transformation = "proportion",
                                   cv_seed = 5466,
                                   nfolds = 10,
                                   queue = "students.q",
                                   dry_lab_subset = NULL){
  
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

  # set up batchtools registry if cross-validation registry does not already exist 
  # and results registry does not already exist
  if(is.null(dry_lab_subset)){
    cv_dir <- paste("registry_paper_glmnet",transformation, sep = "_")
    train_dir <- paste("registry_glmnet_train",transformation,sep = "_")
  } else{
    cv_dir <- paste("registry_paper_glmnet",transformation,
                    dry_lab_subset,sep = "_")
    train_dir <- paste("registry_glmnet_train",transformation,
                       dry_lab_subset,sep = "_")
  }

  if(!(dir.exists(cv_dir)|
                  dir.exists(train_dir))
                 ){ #do not run if results already calculated
    
    if(is.null(dry_lab_subset)){
    glmnet_reg = makeRegistry(file.dir = paste("registry_paper_glmnet",transformation, sep = "_"),
                              packages = c("magrittr",
                                           "glmnet",
                                           "tidyr",
                                           "dplyr",
                                           # "RDS",
                                           "data.table"),
                              source = "all_functions.R")
    } else{
      glmnet_reg = makeRegistry(file.dir = paste("registry_paper_glmnet",
                                                 transformation,
                                                 dry_lab_subset, sep = "_"),
                                packages = c("magrittr",
                                             "glmnet",
                                             "tidyr",
                                             "dplyr",
                                             # "RDS",
                                             "data.table"),
                                source = "all_functions.R")
    
  }
    
    print(glmnet_params)
    
    if(is.null(dry_lab_subset)){
    print(paste("registry_paper_glmnet",transformation, sep = "_"))
    } else{
      print(paste("registry_paper_glmnet",transformation,
                  dry_lab_subset, sep = "_"))
}
    
    
    

    batchMap(fun = do_one_glmnet_cv,
             wet_lab = as.character(glmnet_params[,1]),
             alpha = as.numeric(glmnet_params[,2]),
             class_var = as.character(glmnet_params[,3]),
             level = as.character(glmnet_params[,4]),
             more.args = list(cv_seed = cv_seed,
                              train_seed = train_seed,
                              transformation = transformation,
                              included_dry_labs = included_dry_labs,
                              dry_lab_subset = dry_lab_subset,
                              nfolds = nfolds),
             reg = glmnet_reg)
    
    submitJobs(reg = glmnet_reg,
               resources = list(memory = "25G",
                                queue = queue))
    
    waitForJobs()
    
    return(paste("Submitted jobs for glmnet parameter selection via cross-validation at", 
           transformation, "scale"))
  } else{
    return("glmnet directory already exists")
  }
}

# function to evaluate cross-validated performance of
# xgboost parameters in parallel via batchtools
xgb_cross_validate <- function(included_wet_labs,
                               transformation = "proportion",
                               cv_seed = 5466,
                               nfolds = 10,
                               dry_lab_subset,
                               resources){
  

  
  
  
  #set up batchtools registry
  if(is.null(dry_lab_subset)){
    cv_dir <- paste("registry_paper_xgb",transformation, sep = "_")
    train_dir <- paste("registry_xgb_train",transformation,sep = "_")
  } else{
    cv_dir <- paste("registry_paper_xgb",transformation,
                    dry_lab_subset,sep = "_")
    train_dir <- paste("registry_xgb_train",transformation,
                       dry_lab_subset,sep = "_")
  }
  
  if(!(dir.exists(cv_dir)|
       dir.exists(train_dir))
  ){  #do not run if results already calculated
    xgb_reg = makeRegistry(file.dir = cv_dir,
                           conf.file = "batchtools_multi.conf.R",
                           packages = c("magrittr",
                                        "xgboost",
                                        "tidyr",
                                        "dplyr",
                                        "RDS",
                                        "data.table"),
                           source = "all_functions.R"
    )
    
    #turn off progress bar (improves computation efficiency)
    options(batchtools.progress = FALSE)
    
    
    print(xgb_params)
    
    batchMap(fun = do_one_xgb_cv,
             wet_lab = as.character(xgb_params[,1]),
             eta = as.numeric(xgb_params[,2]),
             max_depth = as.numeric(xgb_params[,3]),
             subsample = as.numeric(xgb_params[,4]),
             colsample_bytree = as.numeric(xgb_params[,5]),
             class_var = as.character(xgb_params[,6]),
             level = as.character(xgb_params[,7]),
             more.args = list(nfolds = nfolds,
                              cv_seed = cv_seed,
                              train_seed = train_seed,
                              transformation = transformation,
                              dry_lab_subset = dry_lab_subset),
             reg = xgb_reg)
    

    
    submitJobs(reg = xgb_reg,
               resources = resources)
    
    waitForJobs()
    return("Done with xgboost parameter selection via cross-validation.")
  } else{
    return("xgboost directory already exists")
  }
}

### function to collect results (returns data.frame with misclassification by lab, classification task)
get_misclass_results <- function(classifier,
                                 transformation,
                                 directory = NULL,
                                 nresults = 112){
  preds <- vector(112,mode = 'list')
  to_read <- paste("registry_",classifier,"_train_",transformation,"/results/",sep = "")
  if(!is.null(directory)){
    to_read <- paste(directory,"/",to_read,sep = "")
  }
  for(i in 1:nresults){
    print(i)
    preds[[i]] <- try(readRDS(paste(to_read,i,".rds",sep="")))
    
  }
  
  misclass <- data.frame(
    predictor_lab = character(0),
    predictee_lab = character(0),
    taxon = character(0),
    class_var = character(0),
    misclass <- numeric(0),
    classifier = character(0),
    within_misclass <- numeric(0),
    equal_misclass = numeric(0)
  )
  
  misclass$classifier %<>% factor(levels = c("Boosted Tree", "Elastic Net"))
  
  if(classifier == "glmnet"){ 
    for(i in 1:nresults){
      
      if(is.list(preds[[i]])){
      misclass_new <- with(preds[[i]],
                           
                           data.frame(
                             predictor_lab = rep(wet_lab,length(included_wet_labs)) %>%
                               as.character(),
                             predictee_lab = glmnet_misclass$wet_lab %>%
                               as.character(),
                             taxon = rep(taxon,length(included_wet_labs)),
                             class_var = rep(class_var,length(included_wet_labs)),
                             misclass = glmnet_misclass$misclass,
                             classifier = ifelse(classifier == "glmnet","Elastic Net","Boosted Tree"),
                             within_misclass = glmnet_misclass$misclass[
                               as.character(glmnet_misclass$wet_lab) == as.character(wet_lab)
                               ] %>% rep(length(glmnet_misclass$misclass)),
                             equal_misclass = glmnet_preds %>%
                               mutate(true_2 = true_labels) %>%
                               group_by(true_labels, wet_lab) %>%
                               summarize(prop_inc = mean(pred_labels != true_2)) %>%
                               group_by(wet_lab) %>%
                               summarize(lab_prop_inc = mean(prop_inc,na.rm = T)) %>%
                               (function(x) x$lab_prop_inc)
      )
      )
      
      misclass_new$classifier %<>%  factor(levels = c("Boosted Tree", "Elastic Net"))
      
      misclass <- rbind(misclass,misclass_new)
      } else{
        warning(paste("Could not load prediction set ", i, sep = "", collapse = ""))
      }
      
      
      
      
    }
  }
  
  if(classifier == "xgb"){
    for(i in 1:nresults){
      if(is.list(preds[[i]])){
      misclass_new <- with(preds[[i]],
                           
                           data.frame(
                             predictor_lab = rep(wet_lab,length(included_wet_labs)) %>%
                               as.character(),
                             predictee_lab = xgb_misclass$wet_lab %>%
                               as.character(),
                             taxon = rep(taxon,length(included_wet_labs)),
                             class_var = rep(class_var,length(included_wet_labs)),
                             misclass = xgb_misclass$misclass,
                             classifier = ifelse(classifier == "glmnet","Elastic Net","Boosted Tree"),
                             within_misclass = xgb_misclass$misclass[
                               as.character(xgb_misclass$wet_lab) == as.character(wet_lab)
                               ] %>% rep(length(xgb_misclass$misclass)),
                             equal_misclass = xgb_preds %>%
                               mutate(true_2 = true_labels) %>%
                               group_by(true_labels, wet_lab) %>%
                               summarize(prop_inc = mean(pred_labels != true_2)) %>%
                               group_by(wet_lab) %>%
                               summarize(lab_prop_inc = mean(prop_inc,na.rm = T)) %>%
                               (function(x) x$lab_prop_inc))
      )
      
      misclass_new$classifier %<>%  factor(levels = c("Boosted Tree", "Elastic Net"))
      
      misclass <- rbind(misclass,misclass_new)
      } else{
        warning(paste("Could not load prediction set ", i, sep = "", collapse = ""))
      }
    }
    
  }
  misclass$taxon %<>%
    as.character %<>%
    factor(levels = c("OTU","species",
                      "genus","family",
                      "order","class",
                      "phylum"))
  misclass$transformation <- transformation
  return(misclass)
}

### calculate confusion matrix for a given classifier
get_confusion <- function(classifier,
                                 transformation,
                                 directory = NULL,
                          specimens = c("sample10",
                          "sample101",
                          "sample102",
                          "sample103",
                          "sample104",
                          "sample11",
                          "sample12",
                          "sample13",
                          "sample14",
                          "sample2",
                          "sample3",
                          "sample4",
                          "sample55",
                          "sample56",
                          "sample57",
                          "sample59",
                          "sample6",
                          "sample61",
                          "sample63",
                          "sample64",
                          "sample8",
                          "sample9"),
                          specimen_types = c("Fecal artificial colony",
                                             "Human",
                                             "Oral artificial colony",
                                             "Robogut")){
  preds <- vector(112,mode = 'list')
  to_read <- paste("registry_",classifier,"_train_",transformation,"/results/",sep = "")
  if(!is.null(directory)){
    to_read <- paste(directory,"/",to_read,sep = "")
  }
  for(i in 1:112){
    print(i)
    preds[[i]] <- try(readRDS(paste(to_read,i,".rds",sep="")))
    
  }
  
  confusion <- vector(112, mode = "list")
  if(classifier == "glmnet"){ 
    for(i in 1:112){
      
      wet_labs <- unique(preds[[i]]$glmnet_preds$wet_lab)
      confusion_matrices <- lapply(wet_labs, function(x) preds[[i]]$glmnet_preds %>% 
                                     (function(y){ y$true_labels %<>% (function(z) factor(z, levels = 0:ifelse(preds[[i]]$class_var == "specimen_type",4,21))); return(y)}) %>%
                                     (function(y){ y$pred_labels %<>% 
                                         (function(z) factor(z,levels = 0:ifelse(preds[[i]]$class_var == "specimen_type",4,21))); return(y)}) %>%
                                    with(table(true_labels,pred_labels)) %>%
                                    (function(y) y[,order(as.numeric(colnames(y)))])) %>%
        (function(x){ names(x) <- wet_labs; return(x)})
      
    
        for(j in 1:length(confusion_matrices)){
          if(preds[[i]]$class_var == "specimen_type"){
            colnames(confusion_matrices[[j]]) <- rownames(confusion_matrices[[j]]) <- specimen_types
          } else{
            colnames(confusion_matrices[[j]]) <- rownames(confusion_matrices[[j]]) <- specimens
        }
      }
      
      confusion[[i]] <- list(predictor_lab = preds[[i]]$wet_lab,
                             class_var = preds[[i]]$class_var,
                             taxon = preds[[i]]$taxon,
                             classifier = classifier,
                             confusion_matrices = confusion_matrices)
                            
    }
  }
  
  if(classifier == "xgb"){
    for(i in 1:112){
      
      wet_labs <- unique(preds[[i]]$xgb_preds$wet_lab)
      confusion_matrices <- lapply(wet_labs, function(x) preds[[i]]$xgb_preds %>% 
                                     filter(wet_lab == x) %>%
                                     (function(y){ y$true_labels %<>% 
                                         (function(z) factor(z, levels = 0:ifelse(preds[[i]]$class_var == "specimen_type",4,21))); return(y)}) %>%
                                     (function(y){ y$pred_labels %<>% 
                                         (function(z) factor(z,levels = 0:ifelse(preds[[i]]$class_var == "specimen_type",4,21))); return(y)}) %>%
                                     with(table(true_labels,pred_labels)) %>%
                                     (function(y) y[,order(as.numeric(colnames(y)))])) %>%
        (function(x){ names(x) <- wet_labs; return(x)})
      
      for(j in 1:length(confusion_matrices)){
        if(preds[[i]]$class_var == "specimen_type"){
          colnames(confusion_matrices[[j]]) <- rownames(confusion_matrices[[j]]) <- specimen_types
        } else{
          colnames(confusion_matrices[[j]]) <- rownames(confusion_matrices[[j]]) <- specimens
        }
      }
      
      confusion[[i]] <- list(predictor_lab = preds[[i]]$wet_lab,
                             class_var = preds[[i]]$class_var,
                             taxon = preds[[i]]$taxon,
                             classifier = classifier,
                             confusion_matrices = confusion_matrices)
     
    }
    
  }
  return(confusion)
}

