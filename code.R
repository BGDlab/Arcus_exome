# manifest file
manifest=read.table('new_manifest_file.tsv',sep='\t',header=T)

# Manuel said 'contrib-dgd', 'contrib-egrp', 'contrib-dgd,contrib-egrp', they all contain exome data

#############################
# check family relationship #
#############################
manifest_exo=manifest[which(manifest$protocol=='Whole Exome Sequencing (WES)'),]
fam=read.csv('family_relationships.csv')

patid_unique=data.frame(patid=unique(manifest_exo$pat_id),
                        count=NA)
for (i in 1:dim(patid_unique)[1]) {
  patid_unique$count[i]=length(which(manifest_exo$pat_id==patid_unique$patid[i]))
}

patid_unique$bioid=manifest_exo$biosample_id[match(patid_unique$patid,manifest_exo$pat_id)]
patid_unique$filename=NA

patid_unique$famid=NA
patid_unique$relation=NA

fam_allid=unique(c(fam$participant_id,fam$relative_id))

overlap=intersect(patid_unique$patid,fam_allid)  
  
remains=setdiff(patid_unique$patid,overlap)

patid_unique$famid[match(remains,patid_unique$patid)]=c(1:length(remains))

# from famid=3728, all are families
tem1=data.frame(sample=unique(fam$participant_id),
                id=c(1:length(unique(fam$participant_id)))+length(remains))

fam$id=tem1$id[match(fam$participant_id,tem1$sample)]

tem2=list()

for (i in 1:dim(tem1)[1]) {
  tem2[[i]]=unique(unname(unlist(fam[which(fam$id==i+length(remains)),1:2])))
}

for (i in 1:dim(patid_unique)[1]) {
  
  if (as.numeric(is.na(patid_unique$famid[i]))==1) {
    
    tem3=c()
    for (j in 1:length(tem2)) {
      tem3[j]=as.numeric(!is.na(match(patid_unique$patid[i],tem2[[j]])))
    }
    
    patid_unique$famid[i]=tem1[which(tem3==1),2]

  }
}

patid_unique$relation=fam$relative_role[match(patid_unique$patid,fam$relative_id)]

# reorder fam_id
tem3=data.frame(famid=sort(unique(patid_unique$famid)),
                famid_reo=c(1:length(unique(patid_unique$famid))))

patid_unique$famid=tem3$famid_reo[match(patid_unique$famid,tem3$famid)]

# Aaron said we need to remove their families
patid_unique=patid_unique[which(is.na(patid_unique$relation)),c('patid','count','bioid','filename')]

# output the sample ID that will be analyzed
name=c()
for (i in 1:dim(patid_unique)[1]) {
  tem4=manifest_exo$filename[which(manifest_exo$pat_id==patid_unique$patid[i])]
  name=c(name,tem4)
}
name=gsub("\\.vcf.gz$", "", name)

write.table(name,file='exome_id.txt', col.names = F, row.names = F, quote = F, sep='')

write.table(patid_unique,file='patid_unique.txt',col.names = T,row.names = F,sep='\t',quote=F)

#########################################################
# summarize number of genetic variations per individual #
#########################################################

# for Polyphen defined 'probably_damaging' genetic mutations
ids=read.table('exome_id.txt',header=F)
ids=ids$V1

num_var_count_v2=c()
for (i in 1:length(ids)) {
  
  aa=read.table(paste('/extract_var_v2/',ids[i],'.txt',sep=''),header=F)
  num_var_count_v2[i]=dim(aa)[1]
  
}

num_var_count_v2=data.frame(id=ids,
                         count=num_var_count_v2)

######################
# check INFO columns #
######################
var_name=c("Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position",
           "Protein_position","Amino_acids","Codons","Existing_variation","ALLELE_NUM","DISTANCE","STRAND","FLAGS","VARIANT_CLASS","SYMBOL_SOURCE","HGNC_ID",
           "CANONICAL","MANE_SELECT","MANE_PLUS_CLINICAL","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","UNIPROT_ISOFORM","GENE_PHENO","SIFT","PolyPhen","DOMAINS",
           "HGVS_OFFSET","AF","gnomADe_AF","gnomADe_AFR_AF","gnomADe_AMR_AF","gnomADe_ASJ_AF","gnomADe_EAS_AF","gnomADe_FIN_AF","gnomADe_NFE_AF","gnomADe_OTH_AF",
           "gnomADe_SAS_AF","gnomADg_AF","gnomADg_AFR_AF","gnomADg_AMI_AF","gnomADg_AMR_AF","gnomADg_ASJ_AF","gnomADg_EAS_AF","gnomADg_FIN_AF","gnomADg_MID_AF",
           "gnomADg_NFE_AF","gnomADg_OTH_AF","gnomADg_SAS_AF","MAX_AF","MAX_AF_POPS","CLIN_SIG","SOMATIC","PHENO","PUBMED","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS",
           "MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS","pLI_gene_value")

aa=read.table('/QB510QZMP6LEE_19697485914516_b69f06dcb4.txt',sep=' ',header=F)
aa=aa[,-9]

count=c()
for (i in 1:dim(aa)[1]) {
  
  aa2=aa$V8[i]
  
  aa3=unlist(strsplit(aa2, split = ";"))[16]
  
  aa4=unlist(strsplit(aa3, split = ","))
  
  tem1=c()
  for (j in 1:length(aa4)) {
    
  aa5=unlist(strsplit(aa4[j], split = "\\|"))
  
  tem1=c(tem1,aa5[37])
  
  }
  
  count[i]=paste(tem1, collapse = ",")
  
}


#######################################################
# summarize number of genetic variations with ClinVar #
#######################################################

ids=read.table('/exome_id.txt',header=F)
ids=ids$V1

num_var_count_polyphen_clinvar=data.frame(id=ids,count=NA)

for (i in 1:length(ids)) {
  
  aa=read.table(paste('/extract_var_v2/',ids[i],'.txt',sep=''),header=F)
  
  tem1=grep('\\|pathogenic\\|',aa$V8) # as indicated by CLN_SIG
  
  tem2=grep('\\|likely_pathogenic\\|',aa$V8) # as indicated by CLN_SIG

  num_var_count_polyphen_clinvar$count[i]=length(unique(c(tem1,tem2)))
  
}

p <- ggplot(num_var_count_polyphen_clinvar, aes(x=count)) +  geom_histogram()

#######################################
# check allele frequency (gnomADe_AF) #
#######################################

# check 'AF', 'gnomADe_AF','gnomADg_AF','MAX_AF', sometimes there is no AF at all

ids=read.table('exome_id.txt',header=F)
ids=ids$V1

num_var_count_polyphen_clinvar=data.frame(id=ids,count=NA,af_count=NA)

for (i in 1:length(ids)) {
  
  aa=read.table(paste('/extract_var_v2/',ids[i],'.txt',sep=''),header=F)
  
  tem1=grep('\\|pathogenic\\|',aa$V8)
  
  tem2=grep('\\|likely_pathogenic\\|',aa$V8)
  
  tem3=unique(c(tem1,tem2))
  
  num_var_count_polyphen_clinvar$count[i]=length(tem3)
  
  if (length(tem3)==0) {
    
    num_var_count_polyphen_clinvar$af_count[i]=0
    
  } else {
    
    tem4=c()
    for (j in 1:length(tem3)){
      
      aa2=aa$V8[tem3[j]]
      
      aa3=unlist(strsplit(aa2, split = ";"))[16]
      
      aa4=unlist(strsplit(aa3, split = ","))
      
      tem5=c()
      tem6=c()
      tem7=c()
      tem8=c()
      for (k in 1:length(aa4)) {
    
      aa5=unlist(strsplit(aa4[k], split = "\\|"))
      
      tem5=c(tem5,aa5[40]) # # as indicated by AF
      tem6=c(tem6,aa5[41]) # as indicated by gnomADe_AF
      tem7=c(tem7,aa5[50]) # as indicated by gnomADg_AF
      tem8=c(tem8,aa5[61]) # as indicated by MAX_AF
      
      }
      
      tem4[j]=max(unique(c(as.numeric(tem5),as.numeric(tem6),as.numeric(tem7),as.numeric(tem8))),na.rm=T)
      
    }
    num_var_count_polyphen_clinvar$af_count[i]=length(which(tem4>0.01))
  }
  
}

num_var_count_polyphen_clinvar$count_clean = num_var_count_polyphen_clinvar$count - num_var_count_polyphen_clinvar$af_count

num_var_count_polyphen_clinvar$patient_id=as.data.frame(do.call(rbind, strsplit(num_var_count_polyphen_clinvar$id, "_")), stringsAsFactors = FALSE)$V1

p <- ggplot(num_var_count_polyphen_clinvar, aes(x=count_clean)) +  geom_histogram()


#############################################
# map pathogenic genetic mutations to genes #
#############################################
  
ids=read.table('exome_id.txt',header=F)
ids=ids$V1

# we already know that there are 3 genetic mutations at most
gm_polyphen_clinvar_af=data.frame(id=ids,gm1=NA,gm2=NA,gm3=NA)

for (i in 1:length(ids)) {
  
  aa=read.table(paste('/extract_var_v2/',ids[i],'.txt',sep=''),header=F)
  
  tem1=grep('\\|pathogenic\\|',aa$V8)
  
  tem2=grep('\\|likely_pathogenic\\|',aa$V8)
  
  tem3=unique(c(tem1,tem2))
  
  
  if (length(tem3)!=0) {

    
    tem4=c()
    for (j in 1:length(tem3)){
      
      aa2=aa$V8[tem3[j]]
      
      aa3=unlist(strsplit(aa2, split = ";"))[16]
      
      aa4=unlist(strsplit(aa3, split = ","))
      
      tem5=c()
      tem6=c()
      tem7=c()
      tem8=c()
      for (k in 1:length(aa4)) {
        
        aa5=unlist(strsplit(aa4[k], split = "\\|"))
        
        tem5=c(tem5,aa5[40]) # # as indicated by AF
        tem6=c(tem6,aa5[41]) # as indicated by gnomADe_AF
        tem7=c(tem7,aa5[50]) # as indicated by gnomADg_AF
        tem8=c(tem8,aa5[61]) # as indicated by MAX_AF
        
      }
      
      tem4[j]=max(unique(c(as.numeric(tem5),as.numeric(tem6),as.numeric(tem7),as.numeric(tem8))),na.rm=T)
      
    }
    
  }
  
    tem3=tem3[which(tem4<0.01)]
    
    if (length(tem3) <= 3) {
      gm_polyphen_clinvar_af[i,2:4] <- c(tem3, rep(NA, 3 - length(tem3)))
    }
    
}

# map genetic mutations to genes
gm_polyphen_clinvar_af_gene=data.frame(id=ids,gene1=NA,gene2=NA,gene3=NA)
gm_polyphen_clinvar_af_ensg=data.frame(id=ids,gene1=NA,gene2=NA,gene3=NA)

for (i in 1:dim(gm_polyphen_clinvar_af_gene)[1]) {
  if (length(which(is.na(gm_polyphen_clinvar_af[i,2:4])))<3) {
    
    aa=read.table(paste('/extract_var_v2/',ids[i],'.txt',sep=''),header=F)
    
    for (j in 1:length(which(!is.na(gm_polyphen_clinvar_af[i,2:4])))) {
    aa2=aa$V8[gm_polyphen_clinvar_af[i,1+j]]
    aa3=unlist(strsplit(aa2, split = ";"))[16]
    aa4=unlist(strsplit(aa3, split = ","))
    
    tem5=c()
    tem6=c()
    for (k in 1:length(aa4)) {
      
      aa5=unlist(strsplit(aa4[k], split = "\\|"))
      
      tem5=c(tem5,aa5[4]) # gene ID
      tem6=c(tem6,aa5[5]) # ENSG ID
      
    }
    
    gm_polyphen_clinvar_af_gene[i,1+j]=unique(tem5[tem5 != ""])[1]
    gm_polyphen_clinvar_af_ensg[i,1+j]=unique(tem6[tem6 != ""])[1]
    
    }
    
    
  }
}

# for each gene, count and plot the number of subjects carrying genetic mutations located in this gene
genes=unique(c(gm_polyphen_clinvar_af_gene[,2],gm_polyphen_clinvar_af_gene[,3],gm_polyphen_clinvar_af_gene[,4]))
genes=genes[-which(is.na(genes))]

num_sub=data.frame(gene_id=genes,
                   num=NA)

for (i in 1:dim(num_sub)[1]) {
  tem1=which(gm_polyphen_clinvar_af_gene==num_sub$gene_id[i], arr.ind = TRUE)
  num_sub$num[i]=length(unique(tem1[,1]))
}

#######################################
# check individuals with imaging data #
#######################################
library(bigrquery)
bq_auth()

project_id <- if (exists('project_id')) project_id else bq_projects()[1]

query <- "SELECT * FROM arcus.brain_mri"

tb <- bq_project_query(project_id, query)

imaging_sql <- bq_table_download(tb)

results$brain_mri_flag=NA

results$brain_mri_flag=imaging_sql$brain_mri_flag[match(results$patient_id,imaging_sql$pat_id)] # 1 means this subject has imaging data

# 1104 individuals (1159 samples) have pathogenic genetic mutations
subid_gene_mut=unique(num_var_count_polyphen_clinvar$patient_id[which(num_var_count_polyphen_clinvar$count_clean>0)])

# among 411 of 1104 individuals have imaging data

subid_gene_mut_imaging=intersect(results$patient_id[which(results$brain_mri_flag==1)] , subid_gene_mut)
write.table(subid_gene_mut_imaging,file='/test/subid_gene_mut_imaging.txt',quote=F,col.names = F,row.names = F,sep='')

