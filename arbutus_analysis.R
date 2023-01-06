#Arbutus Analysis

library(geiger)
library(arbutus)
library(tidyverse)
library(flipR)
library(parallel)

#Load the tree
tree <- read.tree("species_tree/species_tree/species_timetree.nwk")

#Load the data
Orthogroups <- read_delim("orthogroup_statistics/orthogroup_statistics/Orthogroups.csv", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)
#single couple only
single_copy <- read.delim("orthogroup_statistics/orthogroup_statistics/SingleCopyOrthogroups.txt", header = FALSE)

final_orthogroups <- Orthogroups %>% filter(`...1` %in% single_copy$V1) %>% rename(orthogroup = `...1`) %>%
  transmute(orthogroup = orthogroup,
            anolis_carolinensis = str_extract(Anolis_carolinensis, "ENS..............."),
            astyanax_mexicanus = str_extract(Astyanax_mexicanus, "ENS..............."),
            bos_taurus = str_extract(Bos_taurus, "ENS..............."),
            callithrix_jacchus = str_extract(Callithrix_jacchus, "ENS..............."),
            canis_lupus = str_extract(Canis_lupus, "ENS..............."),
            chinchilla_lanigera = str_extract(Chinchilla_lanigera, "ENS..............."),
            danio_rerio = str_extract(Danio_rerio, "ENS..............."),
            gadus_morhua = str_extract(Gadus_morhua, "ENS..............."),
            gallus_gallus = str_extract(Gallus_gallus, "ENS..............."),
            homo_sapiens = str_extract(Homo_sapiens, "ENS............"),
            macaca_mulatta = str_extract(Macaca_mulatta, "ENS..............."),
            monodelphis_domestica = str_extract(Monodelphis_domestica, "ENS..............."),
            mus_musculus = str_extract(Mus_musculus, "ENS..............."),
            oreochromis_niloticus = str_extract(Oreochromis_niloticus, "ENS..............."),
            ornithorhynchus_anatinus = str_extract(Ornithorhynchus_anatinus, "ENS..............."),
            oryzias_latipes = str_extract(Oryzias_latipes, "ENS..............."),
            ovis_aries = str_extract(Ovis_aries, "ENS..............."),
            rattus_norvegicus = str_extract(Rattus_norvegicus, "ENS..............."),
            sus_scrofa = str_extract(Sus_scrofa, "ENS..............."),
            xenopus_tropicalis = str_extract(Xenopus_tropicalis, "ENS..............."),
            oryctolagus_cuniculus = str_extract(Oryctolagus_cuniculus, "ENS..............."))

#Species FPKM data
Anolis_carolinensis <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Anolis_carolinensis.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Anolis_carolinensis") %>% rename(anolis_carolinensis = "Gene")
Astyanax_mexicanus <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Astyanax_mexicanus.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Astyanax_mexicanus") %>% rename(astyanax_mexicanus = "Gene")
Bos_taurus <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Bos_taurus.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Bos_taurus") %>% rename(bos_taurus = "Gene")
Callithrix_jacchus <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Callithrix_jacchus.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Callithrix_jacchus") %>% rename(callithrix_jacchus = "Gene")
Canis_lupus <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Canis_lupus.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Canis_lupus") %>% rename(canis_lupus = "Gene")
Chinchilla_lanigera <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Chinchilla_lanigera.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Chinchilla_lanigera") %>% rename(chinchilla_lanigera = "Gene")
Danio_rerio <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Danio_rerio.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Danio_rerio") %>% rename(danio_rerio = "Gene")
Gadus_morhua <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Gadus_morhua.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Gadus_morhua") %>% rename(gadus_morhua = "Gene")
Gallus_gallus <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Gallus_gallus.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Gallus_gallus") %>% rename(gallus_gallus = "Gene")
Homo_sapiens <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Homo_sapiens.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Homo_sapiens") %>% rename(homo_sapiens = "Gene")
Macaca_mulatta <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Macaca_mulatta.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Macaca_mulatta") %>% rename(macaca_mulatta = "Gene")
Monodelphis_domestica <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Monodelphis_domestica.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Monodelphis_domestica") %>% rename(monodelphis_domestica = "Gene")
Mus_musculus <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Mus_musculus.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Mus_musculus") %>% rename(mus_musculus = "Gene")
Oreochromis_niloticus <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Oreochromis_niloticus.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Oreochromis_niloticus") %>% rename(oreochromis_niloticus = "Gene")
Ornithorhynchus_anatinus <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Ornithorhynchus_anatinus.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Ornithorhynchus_anatinus") %>% rename(ornithorhynchus_anatinus = "Gene")
Oryctolagus_cuniculus <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Oryctolagus_cuniculus.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Oryctolagus_cuniculus") %>% rename(oryctolagus_cuniculus = "Gene")
Oryzias_latipes <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Oryzias_latipes.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Oryzias_latipes") %>% rename(oryzias_latipes = "Gene")
Ovis_aries <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Ovis_aries.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Ovis_aries") %>% rename(ovis_aries = "Gene")
Rattus_norvegicus <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Rattus_norvegicus.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Rattus_norvegicus") %>% rename(rattus_norvegicus = "Gene")
Sus_scrofa <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Sus_scrofa.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Sus_scrofa") %>% rename(sus_scrofa = "Gene")
Xenopus_tropicalis <- read_tsv("transcriptome_amalgamation/transcriptome_amalgamation/sva_log_tmm_fpkm/tissue_mean/Xenopus_tropicalis.tissue.mean.tsv", col_names = c("Gene", "brain", "heart", "kidney", "liver", "ovary", "testis"), skip = 1) %>%
  pivot_longer(cols = -Gene, names_to = "Organ", values_to = "Xenopus_tropicalis") %>% rename(xenopus_tropicalis = "Gene")

total_df <- final_orthogroups %>% full_join(Anolis_carolinensis) %>% full_join(Astyanax_mexicanus) %>% full_join(Bos_taurus) %>% full_join(Callithrix_jacchus) %>%
  full_join(Canis_lupus) %>% full_join(Chinchilla_lanigera) %>% full_join(Danio_rerio) %>% full_join(Gadus_morhua) %>% full_join(Gallus_gallus) %>%
  full_join(Homo_sapiens) %>% full_join(Macaca_mulatta) %>% full_join(Monodelphis_domestica) %>% full_join(Mus_musculus) %>% full_join(Oreochromis_niloticus) %>%
  full_join(Ornithorhynchus_anatinus) %>% full_join(Oryctolagus_cuniculus) %>% full_join(Oryzias_latipes) %>% full_join(Ovis_aries) %>%
  full_join(Rattus_norvegicus) %>% full_join(Sus_scrofa) %>% full_join(Xenopus_tropicalis) %>% select(c(orthogroup, Organ:Xenopus_tropicalis))

rm(Orthogroups, single_copy, Anolis_carolinensis, Astyanax_mexicanus, Bos_taurus, Callithrix_jacchus, Canis_lupus, Chinchilla_lanigera, Danio_rerio, Gadus_morhua, Gallus_gallus, Homo_sapiens, Macaca_mulatta, Mus_musculus, Oreochromis_niloticus, Ornithorhynchus_anatinus, Oryctolagus_cuniculus, Oryzias_latipes, Ovis_aries, Rattus_norvegicus, Sus_scrofa, Xenopus_tropicalis)

#Split data into organs
brain_df <- total_df %>% filter(Organ == "brain") %>% drop_na(orthogroup) %>% select(-Organ)
heart_df <- total_df %>% filter(Organ == "heart") %>% drop_na(orthogroup) %>% select(-Organ)
kidney_df <- total_df %>% filter(Organ == "kidney") %>% drop_na(orthogroup) %>% select(-Organ)
liver_df <- total_df %>% filter(Organ == "liver") %>% drop_na(orthogroup) %>% select(-Organ)
ovary_df <- total_df %>% filter(Organ == "ovary") %>% drop_na(orthogroup) %>% select(-Organ)
testis_df <- total_df %>% filter(Organ == "testis") %>% drop_na(orthogroup) %>% select(-Organ)

#Standard Error function
standard_error <- function(x) sd(x) / sqrt(length(x))

#Get SE for each gene per body part
brain_SE <- brain_df %>% group_by(orthogroup) %>% summarise(orthogroup, SE = standard_error(across(Anolis_carolinensis:Xenopus_tropicalis)))
heart_SE <- heart_df %>% group_by(orthogroup) %>% summarise(orthogroup, SE = standard_error(across(Anolis_carolinensis:Xenopus_tropicalis)))
kidney_SE <- kidney_df %>% group_by(orthogroup) %>% summarise(orthogroup, SE = standard_error(across(Anolis_carolinensis:Xenopus_tropicalis)))
liver_SE <- liver_df %>% group_by(orthogroup) %>% summarise(orthogroup, SE = standard_error(across(Anolis_carolinensis:Xenopus_tropicalis)))
ovary_SE <- ovary_df %>% group_by(orthogroup) %>% summarise(orthogroup, SE = standard_error(across(Anolis_carolinensis:Xenopus_tropicalis)))
testis_SE <- testis_df %>% group_by(orthogroup) %>% summarise(orthogroup, SE = standard_error(across(Anolis_carolinensis:Xenopus_tropicalis)))

#Now need to flip tables and properly format
format_expr_data <- function (avgdat) {
  temp <- avgdat %>% pull(orthogroup)
  avgdat <- avgdat %>% ungroup() %>% select(!orthogroup)
  dat <- flip(avgdat)
  colnames(dat) <- temp
  res <- dat %>% as.matrix()
  res
}

#Running fitcontinuous
runFC <- function ( dat, SE ){
  fitResults <- vector(mode = "list", length = ncol(dat))
  for(j in 1:ncol(dat)){
    tdf <- treedata(tree, dat[,j], sort = TRUE)
    phy <- tdf$phy
    data <- tdf$data
    fitBM <- fitContinuous(phy, data, SE[[2]][[j]], model = "BM")
    fitOU <- fitContinuous(phy, data, SE[[2]][[j]], model = "OU")
    fitEB <- fitContinuous(phy, data, SE[[2]][[j]], model = "EB")
    aic <- c(fitBM$opt[["aic"]], fitOU$opt[["aic"]], fitEB$opt[["aic"]])
    fit <- ifelse(min(aic) == aic[1], list(c(fitBM, model = "BM")), 
                  ifelse(min(aic) == aic[2], list(c(fitOU, model = "OU")), 
                         list(c(fitEB, model = "EB"))))
    fitResults[j] <- fit
  }
  fitResults
}

model_count <- function (fit) {
  ou = 0
  bm = 0
  eb = 0
  for(f in fit){
    vec <- f
    ifelse(vec$model == "OU", ou <- ou + 1, ifelse(vec$model == "BM", bm <- bm + 1, eb <- eb + 1))
  }
  df <- data.frame(OU = ou, BM = bm, EB = eb)
  b <- df %>% pivot_longer(c(OU, BM, EB), names_to = "model")
  b
}

#running arbutus
run_arb <- function (fits){
  arby <- vector("list", length = length(fits))
  count = 1
  for(f in fits){
    class(f) <- "gfit"
    arby[[count]] <- arbutus(f)
    count = count + 1
  }
  arby_df <- map_df(arby, pvalue_arbutus)
  arby_df
}

total_process <- function (dat_list){
  avgdat <- dat_list[[1]]
  part <- dat_list[[2]]
  SE <- dat_list[[3]]
  exp <- format_expr_data(avgdat)
  fit <- runFC(exp, SE)
  fit_name <- paste0("arbutus/fits/fit_", part)
  saveRDS(fit, file = fit_name)
  df <- model_count(fit)
  aic_name <- paste0("arbutus/AIC/AIC_", part, ".png")
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  ggsave(aic_name)
  #result <- run_arb(fit)
  #rds_name <- paste0("arbutus/pvals/pvals_", part)
  #saveRDS(result, file = rds_name)
  #result %>% select(!m.sig) %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
   # ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  #pval_name <- paste0("arbutus/figures/arbutus_", part, ".png")
  #ggsave(pval_name)
}

br_list <- list(brain_df, "brain", brain_SE)
cb_list <- list(heart_df, "heart", heart_SE)
ht_list <- list(kidney_df, "kidney", kidney_SE)
kd_list <- list(liver_df, "liver", liver_SE)
lv_list <- list(ovary_df, "ovary", ovary_SE)
ts_list <- list(testis_df, "testis", testis_SE)
all_list <- list(br_list, cb_list, ht_list, kd_list, lv_list, ts_list)

#lapply(all_list, total_process)
mclapply(all_list, total_process, mc.cores = 6)