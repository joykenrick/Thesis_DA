## LIMMA ## 
library(limma)
library(ggrepel)

limma_per_cancer_batch1 <- function(input_cancer,input_ctrl,cancer){
  
  one_cancer <- filter(input_cancer,Disease==cancer)
  
  filt_combined <- rbind(input_ctrl,one_cancer)
  #filt_combined <- filter(filt_combined,!Assay %in% c("LTA4H", "GDF15", "PMVK", "CTSC"))
  
  meta_combined <- rbind(UCAN_metadata, IGT_metadata)
  
  #cancer_levels <- c('AML','CLL','LYMPH','MYEL','CRC','LUNGC','GLIOM','BRC', 'CVX','ENDC','OVC','PRC','HPC','PAN','MEL')
  
  wide_dat <-
    filt_combined %>% 
    select(DAid,Assay,NPX) %>% 
    spread(Assay,NPX,-1) 
  
  all_data <- merge(wide_dat,meta_combined, by = 'DAid')
  all_data <-
    all_data %>% 
    mutate(GROUP = ifelse(GROUP == cancer, paste("1_", cancer, sep = ""), paste("0_", GROUP, sep = "")))
  
  design<-model.matrix(~0 + all_data$Age + all_data$Sex + all_data$GROUP)
  colnames(design) <- c('Age', "Female","Male",as.character(cancer))
  
  fit <- all_data %>% 
    select(-DAid,-Sex,-Age,-BMI,-Stage,-Grade,-GROUP) %>%
    t() %>% 
    lmFit(design=design)
  
  
  contrast <- makeContrasts(as.character(cancer), levels=design)
  
  # apply contrast
  contrast_fit<-contrasts.fit(fit, contrast)
  
  # apply empirical Bayes smoothing to the SE
  ebays_fit<-eBayes(contrast_fit)
  
  # summary
  print(summary(decideTests(ebays_fit)))
  # extract DE results
  DE_results<-topTable(ebays_fit, n=ncol(all_data), adjust.method="fdr", confint=TRUE)
  
  DE_results <- DE_results %>% 
    mutate(
      Expression = case_when(logFC >= log(2)& -log10(adj.P.Val) >30 ~ "Up-regulated",
                             logFC <= -log(2) & -log10(adj.P.Val) >30 ~ "Down-regulated",
                             TRUE ~ "Unchanged"))
  top <- 15
  top_up <- DE_results %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(desc(-log10(adj.P.Val)), desc(abs(logFC))) %>% 
    head(top)
  top_down <- DE_results %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(desc(-log10(adj.P.Val)), desc(abs(logFC))) %>% 
    head(top)
  
  top_genes <- rbind(top_up,top_down)
  
  #arrange(-log10(adj.P.Val)))
  #slice(1:top)
  
  volcano_plot <- DE_results %>% 
    ggplot(aes(logFC,-log10(adj.P.Val)))+
    geom_point() +
    theme_bw() +
    ggtitle(paste(cancer))
  
  volcano_plot <- volcano_plot +
    geom_label_repel(data = top_up,
                     mapping = aes(logFC, -log10(adj.P.Val), label = rownames(top_up)),
                     size = 2)
  
  return(list(DE_results = DE_results,
              ebays_fit = ebays_fit,
              top_genes = top_genes,
              top_up = top_up,
              top_down = top_down,
              volcano_plot = volcano_plot))
}

limma_per_cancer_female <- function(input_cancer,input_ctrl,cancer){
  
  one_cancer <- filter(input_cancer,Disease==cancer)
  #  one_cancer <- filter(one_cancer, Age > 0)
  input_ctrl <- filter(input_ctrl,Sex == 'FEMALE')
  filt_combined <- rbind(input_ctrl,one_cancer)
  filt_combined <- filter(filt_combined,Age>0)
  
  #filt_combined <- filter(filt_combined,!Assay %in% c("LTA4H", "GDF15", "PMVK", "CTSC"))
  
  #meta_combined <- rbind(UCAN_metadata, IGT_metadata)
  
  #cancer_levels <- c('AML','CLL','LYMPH','MYEL','CRC','LUNGC','GLIOM','BRC', 'CVX','ENDC','OVC','PRC','HPC','PAN','MEL')
  
  wide_dat <-
    filt_combined %>% 
    select(DAid,Assay,NPX,Disease,Age,Sex) %>% 
    spread(Assay,NPX,-1) 
  
  #all_data <- merge(wide_dat,filt_combined, by = 'DAid')
  all_data <- wide_dat
  all_data <-
    all_data %>% 
    mutate(Disease = ifelse(Disease == cancer, paste("1_", cancer, sep = ""), paste("0_", Disease, sep = "")))
  
  design<-model.matrix(~0 + all_data$Age + all_data$Disease)
  colnames(design) <- c('Age', "Healthy",as.character(cancer))
  
  fit <- all_data %>% 
    select(-DAid,-Sex,-Age,-Disease) %>%
    t() %>% 
    lmFit(design=design)
  
  
  contrast <- makeContrasts(as.character(cancer), levels=design)
  
  # apply contrast
  contrast_fit<-contrasts.fit(fit, contrast)
  
  # apply empirical Bayes smoothing to the SE
  ebays_fit<-eBayes(contrast_fit)
  
  # summary
  print(summary(decideTests(ebays_fit)))
  # extract DE results
  DE_results<-topTable(ebays_fit, n=ncol(all_data), adjust.method="fdr", confint=TRUE)
  
  DE_results <- DE_results %>% 
    mutate(
      Expression = case_when(logFC > 1 & -log10(adj.P.Val) >1.3 ~ "Up-regulated",
                             logFC < 1 & -log10(adj.P.Val) >1.3 ~ "Down-regulated",
                             TRUE ~ "Unchanged"))
  top <- 10
  top_up <- DE_results %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(desc(-log10(adj.P.Val)), desc(abs(logFC))) %>% 
    head(top)
  top_down <- DE_results %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(desc(-log10(adj.P.Val)), desc(abs(logFC))) %>% 
    head(top)
  
  top_genes <- rbind(top_up,top_down)
  
  #arrange(-log10(adj.P.Val)))
  #slice(1:top)
  
  volcano_plot <- DE_results %>% 
    ggplot(aes(logFC,-log10(adj.P.Val)))+
    geom_point() +
    theme_bw() +
    ggtitle(paste(cancer))
  
  volcano_plot <- volcano_plot +
    geom_label_repel(data = top_up,
                     mapping = aes(logFC, -log10(adj.P.Val), label = rownames(top_up)),
                     size = 2)
  
  return(list(DE_results = DE_results,
              ebays_fit = ebays_fit,
              top_genes = top_genes,
              top_up = top_up,
              top_down = top_down,
              volcano_plot = volcano_plot))
}
limma_per_cancer_male <- function(input_cancer,input_ctrl,cancer){
  
  one_cancer <- filter(input_cancer,Disease==cancer)
  #  one_cancer <- filter(one_cancer, Age > 0)
  input_ctrl <- filter(input_ctrl,Sex == 'MALE')
  filt_combined <- rbind(input_ctrl,one_cancer)
  filt_combined <- filter(filt_combined,Age>0)
  
  #filt_combined <- filter(filt_combined,!Assay %in% c("LTA4H", "GDF15", "PMVK", "CTSC"))
  
  #meta_combined <- rbind(UCAN_metadata, IGT_metadata)
  
  #cancer_levels <- c('AML','CLL','LYMPH','MYEL','CRC','LUNGC','GLIOM','BRC', 'CVX','ENDC','OVC','PRC','HPC','PAN','MEL')
  
  wide_dat <-
    filt_combined %>% 
    select(DAid,Assay,NPX,Disease,Age,Sex) %>% 
    spread(Assay,NPX,-1) 
  
  #all_data <- merge(wide_dat,filt_combined, by = 'DAid')
  all_data <- wide_dat
  all_data <-
    all_data %>% 
    mutate(Disease = ifelse(Disease == cancer, paste("1_", cancer, sep = ""), paste("0_", Disease, sep = "")))
  
  design<-model.matrix(~0 + all_data$Age + all_data$Disease)
  colnames(design) <- c('Age', "Healthy",as.character(cancer))
  
  fit <- all_data %>% 
    select(-DAid,-Sex,-Age,-Disease) %>%
    t() %>% 
    lmFit(design=design)
  
  
  contrast <- makeContrasts(as.character(cancer), levels=design)
  
  # apply contrast
  contrast_fit<-contrasts.fit(fit, contrast)
  
  # apply empirical Bayes smoothing to the SE
  ebays_fit<-eBayes(contrast_fit)
  
  # summary
  print(summary(decideTests(ebays_fit)))
  # extract DE results
  DE_results<-topTable(ebays_fit, n=ncol(all_data), adjust.method="fdr", confint=TRUE)
  
  DE_results <- DE_results %>% 
    mutate(
      Expression = case_when(logFC > 1 & -log10(adj.P.Val) >1.3 ~ "Up-regulated",
                             logFC < 1 & -log10(adj.P.Val) >1.3 ~ "Down-regulated",
                             TRUE ~ "Unchanged"))
  top <- 10
  top_up <- DE_results %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(desc(-log10(adj.P.Val)), desc(abs(logFC))) %>% 
    head(top)
  top_down <- DE_results %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(desc(-log10(adj.P.Val)), desc(abs(logFC))) %>% 
    head(top)
  
  top_genes <- rbind(top_up,top_down)
  
  #arrange(-log10(adj.P.Val)))
  #slice(1:top)
  
  volcano_plot <- DE_results %>% 
    ggplot(aes(logFC,-log10(adj.P.Val)))+
    geom_point() +
    theme_bw() +
    ggtitle(paste(cancer))
  
  volcano_plot <- volcano_plot +
    geom_label_repel(data = top_up,
                     mapping = aes(logFC, -log10(adj.P.Val), label = rownames(top_up)),
                     size = 2)
  
  return(list(DE_results = DE_results,
              ebays_fit = ebays_fit,
              top_genes = top_genes,
              top_up = top_up,
              top_down = top_down,
              volcano_plot = volcano_plot))
}


cancer = 'CVX'
cvx_limma_2 <- limma_per_cancer_female(filt_cancer,filt_IGT,cancer)
cancer = 'AML'
aml_limma <- limma_per_cancer_batch3(filt_cancer,filt_IGT,cancer)
cancer = 'CLL'
cll_limma <- limma_per_cancer_batch3(filt_cancer,filt_IGT,cancer)
cancer = 'LYMPH'
lymph_limma <- limma_per_cancer_batch3(filt_cancer,filt_IGT,cancer)

cancer = 'MYEL'
myel_limma <- limma_per_cancer_batch3(filt_cancer,filt_IGT,cancer)
cancer = 'CRC'
crc_limma <- limma_per_cancer_batch3(filt_cancer,filt_IGT,cancer)
cancer = 'LUNGC'
lungc_limma <- limma_per_cancer_batch3(filt_cancer,filt_IGT,cancer)
cancer = 'GLIOM'
gliom_limma <- limma_per_cancer_batch3(filt_cancer,filt_IGT,cancer)

cancer = 'CVX'
cvx_limma <- limma_per_cancer_female(filt_cancer,filt_IGT,cancer)
cancer = 'BRC'
brc_limma <- limma_per_cancer_female(filt_cancer,filt_IGT,cancer)
cancer = 'ENDC'
endc_limma <- limma_per_cancer_female(filt_cancer,filt_IGT,cancer)
cancer = 'OVC'
ovc_limma <- limma_per_cancer_female(filt_cancer,filt_IGT,cancer)

cancer = 'PRC'
prc_limma <- limma_per_cancer_male(filt_cancer,filt_IGT,cancer)



limma_per_cancer_batch3 <- function(input_cancer,input_ctrl,cancer){
  
  one_cancer <- filter(input_cancer,Disease==cancer)
#  one_cancer <- filter(one_cancer, Age > 0)
  filt_combined <- rbind(input_ctrl,one_cancer)
  filt_combined <- filter(filt_combined,Age>0)
  
    #filt_combined <- filter(filt_combined,!Assay %in% c("LTA4H", "GDF15", "PMVK", "CTSC"))
  
  #meta_combined <- rbind(UCAN_metadata, IGT_metadata)
  
  #cancer_levels <- c('AML','CLL','LYMPH','MYEL','CRC','LUNGC','GLIOM','BRC', 'CVX','ENDC','OVC','PRC','HPC','PAN','MEL')
  
  wide_dat <-
    filt_combined %>% 
    select(DAid,Assay,NPX,Disease,Age,Sex) %>% 
    spread(Assay,NPX,-1) 
  
  #all_data <- merge(wide_dat,filt_combined, by = 'DAid')
  all_data <- wide_dat
  all_data <-
    all_data %>% 
    mutate(Disease = ifelse(Disease == cancer, paste("1_", cancer, sep = ""), paste("0_", Disease, sep = "")))
  
  design<-model.matrix(~0 + all_data$Age + all_data$Sex + all_data$Disease)
  colnames(design) <- c('Age', "Female","Male",as.character(cancer))
  
  fit <- all_data %>% 
    select(-DAid,-Sex,-Age,-Disease) %>%
    t() %>% 
    lmFit(design=design)
  
  
  contrast <- makeContrasts(as.character(cancer), levels=design)
  
  # apply contrast
  contrast_fit<-contrasts.fit(fit, contrast)
  
  # apply empirical Bayes smoothing to the SE
  ebays_fit<-eBayes(contrast_fit)
  
  # summary
  print(summary(decideTests(ebays_fit)))
  # extract DE results
  DE_results<-topTable(ebays_fit, n=ncol(all_data), adjust.method="fdr", confint=TRUE)
  
  DE_results <- DE_results %>% 
    mutate(
      Expression = case_when(logFC > 1 & -log10(adj.P.Val) >1.3 ~ "Up-regulated",
                             logFC < 1 & -log10(adj.P.Val) >1.3 ~ "Down-regulated",
                             TRUE ~ "Unchanged"))
  sig_genes <- DE_results %>% 
    filter(Expression != 'Unchanged') %>% 
    arrange(desc(-log10(adj.P.Val)), desc(abs(logFC)))
  sig_genes <- mutate(sig_genes,term = rownames(sig_genes))
  top <- 10
  top_up <- DE_results %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(desc(-log10(adj.P.Val)), desc(abs(logFC))) %>% 
    head(top)
  top_down <- DE_results %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(desc(-log10(adj.P.Val)), desc(abs(logFC))) %>% 
    head(top)
  
  top_genes <- rbind(top_up,top_down)
  
  #arrange(-log10(adj.P.Val)))
  #slice(1:top)
  
  volcano_plot <- DE_results %>% 
    ggplot(aes(logFC,-log10(adj.P.Val)))+
    geom_point() +
    theme_bw() +
    ggtitle(paste(cancer))
  
  volcano_plot <- volcano_plot +
    geom_label_repel(data = top_up,
                     mapping = aes(logFC, -log10(adj.P.Val), label = rownames(top_up)),
                     size = 2)
  
  return(list(DE_results = DE_results,
              sig_genes = sig_genes,
              ebays_fit = ebays_fit,
              top_genes = top_genes,
              top_up = top_up,
              top_down = top_down,
              volcano_plot = volcano_plot))
}
cancer = 'HCC'
HCC_limma <- limma_per_cancer_batch3(filt_cancer,filt_IGT,cancer)
cancer = 'PAN'
PAN_limma <- limma_per_cancer_batch3(filt_cancer,filt_IGT,cancer)
cancer = 'MEL'
MEL_limma <- limma_per_cancer_batch3(filt_cancer,filt_IGT,cancer)

# I filtered for just early stage lung then ran for one cancer at a time # 
cancer = 'LUNGC'
lung_earlystage <- limma_per_cancer_batch3(early_stage_filt_cancer,filt_IGT_Stage,'LUNGC')

cancer = 'CRC'
crc_earlystage <- limma_per_cancer_batch3(early_stage_filt_cancer,filt_IGT_Stage,'CRC')


limma_per_stage_per_cancer <- function(input_cancer,input_ctrl,cancer){
  
  one_cancer <- filter(filt_cancer_withStage,Disease==cancer)
  #  one_cancer <- filter(one_cancer, Age > 0)
  filt_combined <- rbind(filt_IGT_Stage,one_cancer)
  filt_combined <- filter(filt_combined,Age>0)
  filt_combined <- filter(filt_combined,Stage!='Late')
  #filt_combined <- filter(filt_combined,!Assay %in% c("LTA4H", "GDF15", "PMVK", "CTSC"))
  
  #meta_combined <- rbind(UCAN_metadata, IGT_metadata)
  
  #cancer_levels <- c('AML','CLL','LYMPH','MYEL','CRC','LUNGC','GLIOM','BRC', 'CVX','ENDC','OVC','PRC','HPC','PAN','MEL')
  
  wide_dat <-
    filt_combined %>% 
    select(DAid,Assay,NPX,Disease,Age,Sex,Stage) %>% 
    spread(Assay,NPX,-1) 
  
  #all_data <- merge(wide_dat,filt_combined, by = 'DAid')
  all_data <- wide_dat
  all_data <-
    all_data %>% 
    #Dropping NAs but I might want to restore later 
    drop_na() %>% 
    mutate(Disease = ifelse(Disease == cancer, paste("1_", cancer, sep = ""), paste("0_", Disease, sep = "")))
  
  design<-model.matrix(~0 + all_data$Age + all_data$Sex + all_data$Disease + all_data$Stage)
  colnames(design) <- c('Age', "Female","Male",as.character(cancer),'Healthy')
  
  fit <- all_data %>% 
    select(-DAid,-Sex,-Age,-Disease,-Stage) %>%
    t() %>% 
    lmFit(design=design)
  
  
  contrast <- makeContrasts(as.character(cancer), levels=design)
  
  # apply contrast
  contrast_fit<-contrasts.fit(fit, contrast)
  
  # apply empirical Bayes smoothing to the SE
  ebays_fit<-eBayes(contrast_fit)
  
  # summary
  print(summary(decideTests(ebays_fit)))
  # extract DE results
  DE_results<-topTable(ebays_fit, n=ncol(all_data), adjust.method="fdr", confint=TRUE)
  
  DE_results <- DE_results %>% 
    mutate(
      Expression = case_when(logFC > 1 & -log10(adj.P.Val) >1.3 ~ "Up-regulated",
                             logFC < 1 & -log10(adj.P.Val) >1.3 ~ "Down-regulated",
                             TRUE ~ "Unchanged"))
  top <- 10
  top_up <- DE_results %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(desc(-log10(adj.P.Val)), desc(abs(logFC))) %>% 
    head(top)
  top_down <- DE_results %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(desc(-log10(adj.P.Val)), desc(abs(logFC))) %>% 
    head(top)
  
  top_genes <- rbind(top_up,top_down)
  
  #arrange(-log10(adj.P.Val)))
  #slice(1:top)
  
  volcano_plot <- DE_results %>% 
    ggplot(aes(logFC,-log10(adj.P.Val)))+
    geom_point() +
    theme_bw() +
    ggtitle(paste(cancer))
  
  volcano_plot <- volcano_plot +
    geom_label_repel(data = top_up,
                     mapping = aes(logFC, -log10(adj.P.Val), label = rownames(top_up)),
                     size = 2)
  
  return(list(DE_results = DE_results,
              ebays_fit = ebays_fit,
              top_genes = top_genes,
              top_up = top_up,
              top_down = top_down,
              volcano_plot = volcano_plot))
}
