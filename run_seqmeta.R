##########################################################################################
################# check if the required packages have been installed #####################
##########################################################################################

is.installed <- function(requirePackage){ 
  is.element(requirePackage, installed.packages()[, 1])
}

if(!is.installed("rbgen")) {
  install.packages("http://www.well.ox.ac.uk/~gav/resources/rbgen_v1.1.5.tgz", repos = NULL, type = "source")
}

if(!is.installed("seqMeta")) {
  install.packages("seqMeta")
}

library(seqMeta, quietly=TRUE)
library(rbgen,quietly=TRUE)
library(data.table,quietly=TRUE)
library(tictoc)
##########################################################################################
###################################### argument parsing ##################################
##########################################################################################
bgen_file <- ""
snp_info_file <- ""
phenotype_file <- ""
out_prefix <- "seqmate_out"
region <- ""
region_given <- FALSE
region_chr <- ""
region_begin <- -1
region_end <- -1
region_gene <- ""
target_gene <- ""
current_model <- ""
model_format <- ""
current_study <- ""
kb <- 0

theHelpMessge =
  paste(" The program required packages are \"seqMeta\" and \"rbgen\".\n",
        " The usage is: \n",
        " \"Rscript --vanilla run_seqmeta_dev.R --bgen <bgen> --snp <snpINFO_File> --pheno <phenotype_File> --region <chr:begin-end> [--out <prefix>]\n",
        " Arguments: \n",
        "	--bgen		bgen file.\n",
        "	--snp		SNP INFO file.\n",
	"	--pheno		phenotype file.\n",
	"	--out		(optional), output prefix includeing the paht. Default is seqmate_out.\n",
	"	--region	(optional), region. Format: chr:begin-end.\n",
	"       --gene		(optional),gene name. Either --region or --gene is required, and they are mutually exclusive\n",
	"       --window	(optional),window size,in KB,all variants with physical position no more than half the specified kb distance (decimal permitted) from the named gene  are loaded as well.\n",
	"       --model         M1 has PC1-3; M2 has sex, aaoaae(age of last visit or age at onset), PC1-3; M3 has APOEdose, PC1-3; M4 has sex, aaoaae,APOEdose, PC1-3.\n",
        "	-h or --help,	print the help messgae.\n", sep = "")
args <- commandArgs(TRUE)
if(length(args) < 1) {
  args <- c("--help")
}
## Help section
if("--help" %in% args | "-h" %in% args) {
  cat(theHelpMessge)
  q(save="no")
}

### optional ###
if("--out" %in% args){
  argIndex <- which(args == "--out")
  out_prefix <- args[argIndex + 1]
}

### required ###
if("--bgen" %in% args){
  argIndex <- which(args == "--bgen")
  bgen_file <- args[argIndex + 1]
}
if("--snp" %in% args){
  argIndex <- which(args == "--snp")
  snp_info_file <- args[argIndex + 1]
}
if("--pheno" %in% args){
  argIndex <- which(args == "--pheno")
  phenotype_file <- args[argIndex + 1]
}
if("--region" %in% args){
  argIndex <- which(args == "--region")
  region <- args[argIndex + 1]
}

if("--gene" %in% args){
  argIndex <- which(args == "--gene")
  target_gene <- args[argIndex + 1]
}

if("--window" %in% args){
  argIndex <- which(args == "--window")
  kb <- args[argIndex + 1]
}

if("--model" %in% args){
  argIndex <- which(args == "--model")
  current_model<- args[argIndex + 1]
}

if("--study" %in% args){
  argIndex <- which(args == "--study")
  current_study<- args[argIndex + 1]
}
out_prefix
snp_info_file
phenotype_file
region
target_gene
kb
out_prefix


if(bgen_file == "" | snp_info_file == "" | phenotype_file == "" | (region == "" & target_gene == "") | (region != "" & target_gene != "")){
  cat(theHelpMessge)
  q(save="no")
}

if(target_gene != "" & region != ""){
  cat("Please use gene name OR gene region")
  q(save="no")
}


if(region != ""){
  region_chr <- strsplit(region, "[\\:-]+")[[1]][1]
  region_begin <- as.numeric(strsplit(region, "[\\:-]+")[[1]][2])
  region_end <- as.numeric(strsplit(region, "[\\:-]+")[[1]][3])
  if(!is.null(region_chr) & !is.na(region_begin) & !is.na(region_end) & region_begin < region_end & region_end != 0){
    region_given <- TRUE
    if(kb !=0){
      region_begin = region_begin - as.numeric(kb)*1000
      region_end = region_end + as.numeric(kb)*1000
    }
  }
}


message(paste("The input bgen_file is: ", bgen_file, sep = ''))
message(paste("The input SNP INFO file is: ", snp_info_file, sep = ''))
message(paste("The input phenotype file is: ", phenotype_file, sep = ''))
message(paste("The output prefix is: ", out_prefix, sep = ""))
if(region_given){
  message(paste("The region is: ", region_chr, ":" , region_begin, "-", region_end, sep = ''))
} else{
  message(paste("The target gene is: ", target_gene, sep = ''))
}
if (current_model==""){
  cat("Please enter the model name")
  q(save="no")
} else if (current_model=="M1"){
  model_format="status~pc1+pc2+pc3"
  message("The model is adjusted by PC1,PC2,PC3")
} else if (current_model=="M2"){
  model_format="status~aaoaae+sex+pc1+pc2+pc3"
  message("The model is adjusted by AGE,SEX,PC1,PC2,PC3")
} else if (current_model=="M3"){
  model_format="status~apoe4dose+pc1+pc2+pc3"
  message("The model is adjusted by APOE_DOSE,PC1,PC2,PC3")
} else if (current_model=="M4"){
  model_format="status~aaoaae+sex+apoe4dose+pc1+pc2+pc3"
  message("The model is adjusted by AGE,SEX,APOE_DOSE,PC1,PC2,PC3")
}
##########################################################################################
######################################## load snpinfo ####################################
##########################################################################################
date()
tic("Completed")
tic("loaded snpinfo")
message("1/6 Load snpinfo")
snpinfo <- read.table(snp_info_file, header=T, nrows = 5)

if (ncol(snpinfo) > 4){ # extended
   if (region_given){
      snpinfo <- fread(file=snp_info_file, header=T, nThread=2, colClasses=c("character", "character", "numeric", "NULL", "NULL", "character"))
      #snpinfo <- fread(file=snp_info_file, header=T, nThread=2)
   }else{
      snpinfo <- fread(cmd=paste("zgrep", "-w", target_gene, snp_info_file),
                       header=F, nThread=2,
                       colClasses=c("character", "character", "numeric", "NULL", "NULL", "character"), 
                      )
      #snpinfo <- fread(cmd=paste("zgrep", "-w", target_gene, snp_info_file),header=F,nThread=2)
   }
   colnames(snpinfo) <- c("Name", "chr", "position", "gene")
   #colnames(snpinfo) <- c("Name", "gene")
   toc()
}else{
    if (target_gene!=""){
       snpinfo <- fread(cmd=paste("zgrep", "-w", target_gene, snp_info_file), header=F, nThread=2, colClasses=c("character", "character", "character", "character"), col.names=c("Name","chromosome", "position", "gene"))
      # snpinfo <- fread(cmd=paste("zgrep", "-w", target_gene, snp_info_file), header=F, nThread=2, colClasses=c("character", "character"))
    }else{
       snpinfo <- fread(snp_info_file, header=T, colClasses=c("character", "character", "character", "character"), col.names=c("Name","chromosome", "position", "gene"))
    }
    #colnames(snpinfo) <- c("Name", "gene")
    toc()
    tic("formatting")
    #snpinfo$chr <- sapply(strsplit(snpinfo$Name,"\\:"),"[",1)
    #snpinfo$position <- as.numeric(sapply(strsplit(snpinfo$Name,"\\:"),"[",2))
    snpinfo <- setDT(snpinfo)[, c("chr", "position") := tstrsplit(Name, ":", type.convert = T, fixed = T, keep=c(1,2))]
    toc()
}

#snpinfo <- read.table(snp_info_file, header=T)
#colnames(snpinfo) <- c("Name","gene")
#snpinfo$Name = as.character(snpinfo$Name)
#snpinfo$gene = as.character(snpinfo$gene)

# Only keep SNPs that are in the given region
#snpinfo$chr <- sapply(strsplit(snpinfo$Name,"\\:"),"[",1)
#snpinfo$position <- as.numeric(sapply(strsplit(snpinfo$Name,"\\:"),"[",2))
date()

# Only keep SNPs that are in the given region
if(region_given){ 
  tic("region filter")
  target_snpinfo <- subset(snpinfo, chr==region_chr & position >= region_begin & position <= region_end)
  region_gene <- names(which.max(table(snpinfo$gene)))
  region_gene

} else if(target_gene !=""){
  tic("target gene filter")
  target_snpinfo <- snpinfo[snpinfo$gene==target_gene,]
  dim(target_snpinfo)
  if (dim(target_snpinfo)[1] == 0){
    print(paste("ERROR: ", target_gene, " is not found in ", snp_info_file, ".", sep=''))
    q(save="no")
  } else {
    region_given <- TRUE
    region_chr <- target_snpinfo$chr[1]
    region_begin <- min(target_snpinfo$position)
    region_end <- max(target_snpinfo$position)
    region=paste0(region_chr,":",region_begin,"-",region_end)
    region_gene = target_gene
    target_snpinfo <- subset(snpinfo, chr==region_chr & position >= region_begin & position <= region_end)
    if(kb !=0){
      region_begin = region_begin - as.numeric(kb)*1000
      region_end = region_end + as.numeric(kb)*1000
      target_snpinfo <- subset(snpinfo, chr==region_chr & position >= region_begin & position <= region_end)
    }
    region=paste0(region_chr,":",region_begin,"-",region_end)
    print(paste("The region for ", target_gene, " is ", region_chr, ":", region_begin, "-", region_end, sep=''))
    
  }
  dim(target_snpinfo)
}
toc()
tic("2/6 Load phenotype")
phenotype <- read.table(phenotype_file, header=T,sep="\t")
phenotype$status <- phenotype$status - 1 # [1,2] --> [0,1]


##########################################################################################
################################### load dosage from bgen ################################
##########################################################################################
toc()
if (region_given){
  tic("3/6 Load bgen - region")
  date()
  bgen_db = bgen.load(bgen_file, data.frame(chromosome = region_chr, start = region_begin, end = region_end),index.filename = sprintf( "%s.bgi", bgen_file))
  date()
  genotype=t((bgen_db$data[,,'g=1']*1+bgen_db$data[,,'g=2']*2))

  if(length(bgen_db$variants$rsid) <= 1){
  print(paste("ERROR: Need at least two variants. ", length(bgen_db$variants$rsid), " variant is detected.", sep=''))
  q(save="no")
  }
  
  colnames(genotype) <- as.character(bgen_db$variants$rsid)
  rownames(genotype) <- NULL
  genotype = as.matrix(genotype)
  mode(genotype) <- "numeric"
  toc()

  tic("4/6 Remove all subjects with non-[0,1] status")
  good_status_idx <- which( phenotype$status >=0 & phenotype$status <=1 )
  nsubjects_phenotype <- nrow(phenotype) # number of subject in the input
  nsubjects_genotype <- nrow(genotype)
  phenotype <- phenotype[good_status_idx,]
  genotype <- genotype[good_status_idx,]
  nsubjects_good <- nrow(phenotype) 
  
  dim(target_snpinfo)
  toc()
  tic("5/6 Run prepScores")
  message(nrow(target_snpinfo))

  message("running prepScores")
  date()
  #prepScores_results <- prepScores(Z=genotype, status~pc1+pc2+pc3, family=binomial(), SNPInfo=target_snpinfo, data=phenotype)
  prepScores_results <- prepScores(Z=genotype, model_format, family=binomial(), SNPInfo=target_snpinfo, data=phenotype)
  date()
  message("running singlesnpMeta")
  singleSnp_results <- singlesnpMeta(prepScores_results, SNPInfo=target_snpinfo, studyBetas=TRUE)
  
  if (region_given) {
    segMeta_output_RData <- paste0(out_prefix,"/",current_study,"/",current_model,"/",region_chr,"/RData/")

    segMeta_output_RData

    dir.create(segMeta_output_RData, showWarnings=F, recursive=T, mode = "0755")
    saveRDS(prepScores_results, file=paste0(segMeta_output_RData,region_gene,".RData"))

#     segMeta_output_txt <- paste0(out_prefix,"/",region_chr, "/txt/")
#     dir.create(segMeta_output_txt, showWarnings=F, recursive=T, mode = "0755")
#     capture.output(prepScores_results, file=paste0(segMeta_output_txt,region_gene, ".txt"))

    segMeta_output_single <- paste0(out_prefix,"/",current_study,"/",current_model,"/",region_chr, "/single/")
    dir.create(segMeta_output_single, showWarnings=F, recursive=T, mode = "0755")
    write.table(singleSnp_results, file=paste0(segMeta_output_single,region_gene, ".txt"),quote=F,col.names=T,row.names=F)
  } else {
    segMeta_output_RData <- paste0(out_prefix,"/",current_study,"/",current_model,"/",region_chr,"/RData/")
    dir.create(segMeta_output_RData, showWarnings=F, recursive=T, mode = "0755")
    saveRDS(prepScores_results, file=paste0(segMeta_output_RData,target_gene,".RData"))

#     segMeta_output_txt <- paste0(out_prefix,"/",region_chr, "/txt/")
#     dir.create(segMeta_output_txt, showWarnings=F, recursive=T, mode = "0755")
#     capture.output(prepScores_results, file=paste0(segMeta_output_txt,target_gene, ".txt"))

    segMeta_output_single <- paste0(out_prefix,"/",current_study,"/",current_model,"/",region_chr, "/single/")
    dir.create(segMeta_output_single, showWarnings=F, recursive=T, mode = "0755")
    write.table(singleSnp_results, file=paste0(segMeta_output_single,target_gene, ".txt"),quote=F,col.names=T,row.names=F)
    date()
  }
  toc()

} else if(region_chr !="" &  is.na(region_begin) & is.na(region_end)){
   
  genes=unique(snpinfo$gene)
  for(i in c(1:length(genes)))
    {
      tmp_snpinfo=snpinfo[snpinfo$gene==genes[i],]
      tic("3/6 Load bgen - gene") 
      date()
      bgen_db = bgen.load(bgen_file,rsids=as.vector(snpinfo[snpinfo$gene==genes[i],]$Name))
      toc()
      date()
      genotype=t((bgen_db$data[,,'g=1']*1+bgen_db$data[,,'g=2']*2))

      if(length(bgen_db$variants$rsid) <= 1){
        print(paste("ERROR: Need at least two variants. ", length(bgen_db$variants$rsid), " variant is detected.", sep=''))
        q(save="no")
      }

      colnames(genotype) <- as.character(bgen_db$variants$rsid)
      rownames(genotype) <- NULL
      genotype = as.matrix(genotype)
      mode(genotype) <- "numeric"
      
      tic("4/6 Remove all subjects with non-[0,1] status")
      good_status_idx <- which( phenotype$status >=0 & phenotype$status <=1 )
      nsubjects_phenotype <- nrow(phenotype) # number of subject in the input
      nsubjects_genotype <- nrow(genotype)
      phenotype <- phenotype[good_status_idx,]
      genotype <- genotype[good_status_idx,]
      nsubjects_good <- nrow(phenotype)
      toc()

      tic("5/6 Run prepScores")
      date()
      prepScores_results <- prepScores(Z=genotype, model_format, family=binomial(), SNPInfo=tmp_snpinfo, data=phenotype)
      date()
      singleSnp_results <- singlesnpMeta(prepScores_results, SNPInfo=tmp_snpinfo, studyBetas=TRUE)

      segMeta_output_RData <- paste0(out_prefix,"/",current_model,"/",region_chr,"/RData/")
      dir.create(segMeta_output_RData, showWarnings=F, recursive=T, mode = "0755")
      saveRDS(prepScores_results, file=paste0(segMeta_output_RData,genes[i],".RData"))

#       segMeta_output_txt <- paste0(out_prefix,"/",region_chr, "/txt/")
#       dir.create(segMeta_output_txt, showWarnings=F, recursive=T, mode = "0755")
#       capture.output(prepScores_results, file=paste0(segMeta_output_txt,genes[i], ".txt"))

      segMeta_output_single <- paste0(out_prefix,"/",current_study,"/",current_model,"/",region_chr, "/single/")
      dir.create(segMeta_output_single, showWarnings=F, recursive=T, mode = "0755")
      write.table(singleSnp_results, file=paste0(segMeta_output_single,genes[i], ".txt"),quote=F,col.names=T,row.names=F)
      toc()
   
    }
 }

print(gc())
toc()
