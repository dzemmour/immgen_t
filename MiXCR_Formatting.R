library(readxl)
library(xlsx)

sc <- readRDS("seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering_fullTCRinfo_productiveonly.Rds")

maj_clonotypes <- names(which(sort(table(sc$clonotype_id), decreasing = T)>2))

for(a in 1:length(maj_clonotypes)) {
  sample_maj <- names(which(table(sc$clonotype_id, sc$hash.ID)[maj_clonotypes[a],]>0))
  
  all_c <- c()
  
  for(j in 1:length(sample_maj)){
    masterID <- read_xlsx("MiXCR/MiXCR_commands.xlsx", sheet = sample_maj[j],col_names = F)
    
    cells <- masterID$...1[masterID$...2%in%substr(rownames(sc@meta.data)[which(sc$clonotype_id==maj_clonotypes[a] & sc$hash.ID==sample_maj[j])],1,16)]
    bcs <- masterID$...2[masterID$...2%in%substr(rownames(sc@meta.data)[which(sc$clonotype_id==maj_clonotypes[a] & sc$hash.ID==sample_maj[j])],1,16)]

    for(i in 1:length(cells)){
      
      print(cells[i])
      c <- read.table(paste0("MiXCR/",sample_maj[j],"/clones/",cells[i],"_clones.txt"),sep="\t",header = T)
      c <- c[,c('cloneCount','allVHitsWithScore','allDHitsWithScore','allJHitsWithScore','nSeqCDR3','aaSeqCDR3','refPoints')]
      c <- c[c$cloneCount>4,]
      
      if (nrow(c) == 0) {
        next
      }
      
      c$condition <- sample_maj[j]
      c$barcode <- masterID$...2[masterID$...1==cells[i]]
      c$cellID <- cells[i]
      c$nUMI <- as.vector(sc$nCount_RNA[bcs[i]])
      c <- c[,c('barcode','cellID','nUMI','condition','cloneCount','allVHitsWithScore','allDHitsWithScore','allJHitsWithScore','nSeqCDR3','aaSeqCDR3','refPoints')]
      c[nrow(c)+1,] <- ""
      
      if(length(all_c)>length(c)){
        for(l in 12:length(all_c)){
          c[,colnames(all_c)[l]]<-NA
        }
        all_c<-rbind(all_c,c)
      }else{
        all_c<-rbind(all_c,c)
      }
      
    }

    all_c$SeqMatch<-0
    
    for(s in 1:length(all_c$SeqMatch)){
      if(all_c$nSeqCDR3[s]==""){
        # Checks for empty rows
        all_c$SeqMatch[s] <- 0
      }else{
        all_c$SeqMatch[s] <- as.vector(table(all_c$nSeqCDR3)[which(names(table(all_c$nSeqCDR3))==all_c$nSeqCDR3[s])])
      }
    }

    all_c$vd_nregion <- rep(NA,length(all_c$refPoints))
    all_c$vregion <- rep(NA,length(all_c$refPoints))
    all_c$dj_nregion <- rep(NA,length(all_c$refPoints))
    all_c$d_region <- ""
    all_c$allDHitsWithScore[is.na(all_c$allDHitsWithScore)] <- ""

    for(b in 1:length(all_c$refPoints)){
      if(all_c$barcode[b]!=""){
        refpoints_legend <- c('cdr3begin','v3deletion','vend','dbegin','d5deletion','d3deletion','dend','jbegin','j5deletion','cd3end')
        refpoints_numbers <- strsplit(all_c$refPoints,":")[[b]][-c(1:9)][1:10]
        all_c[b,refpoints_legend] <- as.numeric(refpoints_numbers)
        all_c$vregion[b] <- substr(all_c$nSeqCDR3[b],all_c$cdr3begin[b]+1,all_c$vend[b])
        
        if(as.numeric(all_c$v3deletion[b])<=0){
          if(all_c$jbegin[b]!=all_c$vend[b]){
            if(substr(all_c$allVHitsWithScore[b],1,4)=="TRBV" & all_c$allDHitsWithScore[b]!=""){
              if(all_c$jbegin[b]!=all_c$dend[b]){
                if(all_c$j5deletion[b]>0){
                  all_c$dj_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$dend[b]+1,all_c$jbegin[b]+all_c$j5deletion[b])
                }else{
                  all_c$dj_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$dend[b]+1,all_c$jbegin[b])
                }
              }else{
                if(all_c$d3deletion[b]>0){
                  all_c$d_region[b]<-substr(all_c$nSeqCDR3[b],all_c$dbegin[b]+1,all_c$dend[b]-all_c$d3deletion[b])
                  all_c$dj_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$dend[b]-all_c$d3deletion[b]+1,all_c$jbegin[b])
                }
              }
              if(all_c$vend[b]!=all_c$dbegin[b]){
                all_c$vd_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$vend[b]+1,all_c$dbegin[b])
              }
            }else{
              all_c$vd_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$vend[b]+1,all_c$jbegin[b])
            }
          }else{
            all_c$j_region[b]<-substr(all_c$nSeqCDR3[b],all_c$jbegin[b]+1,all_c$cd3end[b])
          }
        }else{
          all_c$vregion[b]<-substr(all_c$nSeqCDR3[b],as.numeric(all_c$cdr3begin[b])+1,as.numeric(all_c$vend[b])-as.numeric(all_c$v3deletion[b]))
          if(is.na(all_c$dbegin[b])){
            all_c$dbegin[b]<-Inf
          }
          all_c$vd_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$vend[b]-all_c$v3deletion[b]+1,min(all_c$dbegin[b],all_c$jbegin[b]))
        }

        if(substr(all_c$allVHitsWithScore[b],1,4)=="TRBV" & all_c$allDHitsWithScore[b]!=""){
          if(all_c$d5deletion[b]<=0 & all_c$d3deletion[b]<=0){
            all_c$d_region[b]<-substr(all_c$nSeqCDR3[b],all_c$dbegin[b]+1,all_c$dend[b])
          }else{
            if(all_c$d5deletion[b]>0){
              if(all_c$d3deletion[b]>0){
                all_c$d_region[b]<-substr(all_c$nSeqCDR3[b],all_c$dbegin[b]+1+all_c$d5deletion[b],all_c$dend[b]-all_c$d3deletion[b])
                all_c$vd_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$dbegin[b]+1,all_c$dbegin[b]+all_c$d5deletion[b])
              }else{
                if(all_c$dbegin[b]-all_c$vend[b]>0){
                  all_c$vd_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$vend[b]+1,all_c$dbegin[b])
                  if(all_c$d5deletion[b]>0){
                    if(all_c$d3deletion[b]>0){
                      all_c$vd_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$vend[b]+1,all_c$dbegin[b]+all_c$d5deletion[b])
                    }
                    all_c$vd_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$vend[b]+1,all_c$dbegin[b]+all_c$d5deletion[b])
                  }
                }else{
                  all_c$vd_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$dbegin[b]+1,all_c$dbegin[b]+all_c$d5deletion[b])
                  if(all_c$dbegin[b]+all_c$d5deletion[b]>all_c$jbegin[b]){
                    all_c$vd_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$vend[b]+1,all_c$jbegin[b])
                    all_c$dj_nregion[b]<-""
                  }
                }
                all_c$d_region[b]<-substr(all_c$nSeqCDR3[b],all_c$dbegin[b]+1+all_c$d5deletion[b],all_c$dend[b])
              }
            }else{
              if(all_c$d3deletion[b]>0){
                all_c$d_region[b]<-substr(all_c$nSeqCDR3[b],all_c$dbegin[b]+1,all_c$dend[b]-all_c$d3deletion[b])
              }else{
                all_c$d_region[b]<-substr(all_c$nSeqCDR3[b],all_c$dbegin[b],all_c$dend[b]-all_c$d3deletion[b])
              }
            }
          }
          if(all_c$dend[b]!=all_c$jbegin[b]){
            if(all_c$dbegin[b]+all_c$d5deletion[b]<all_c$jbegin[b]){
              if(all_c$j5deletion[b]>0){
                all_c$dj_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$dend[b]+1,all_c$jbegin[b]+all_c$j5deletion[b])
              }else{
                all_c$dj_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$dend[b]+1,all_c$jbegin[b])
              }

            }
            if(all_c$d5deletion[b]<=0){
              all_c$d_region[b]<-substr(all_c$nSeqCDR3[b],all_c$dbegin[b]+1,all_c$dend[b])
            }else{
              all_c$d_region[b]<-substr(all_c$nSeqCDR3[b],all_c$dbegin[b]+all_c$d5deletion[b]+1,all_c$dend[b])
            }
          }
        }


        if(all_c$j5deletion[b]<=0){
          all_c$j_region[b]<-substr(all_c$nSeqCDR3[b],all_c$jbegin[b]+1,all_c$cd3end[b])
        }else{
          if(all_c$j5deletion[b]>0){
            all_c$j_region[b]<-substr(all_c$nSeqCDR3[b],all_c$jbegin[b]+all_c$j5deletion[b]+1,all_c$cd3end[b])
            if((all_c$dj_nregion[b]==""|is.na(all_c$dj_nregion[b]))){
              all_c$dj_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$jbegin[b]+1,all_c$jbegin[b]+all_c$j5deletion[b])
            }
          }
        }
        
        if(all_c$jbegin[b]<all_c$vend[b]){
          all_c[b,'dj_nregion']<-""
          if(all_c$v3deletion[b]>0){
            all_c$vregion[b]<-substr(all_c$nSeqCDR3[b],all_c$cdr3begin[b]+1,all_c$vend[b]-all_c$v3deletion[b])
            all_c$vd_nregion[b]<-substr(all_c$nSeqCDR3[b],all_c$vend[b]-all_c$v3deletion[b]+1,all_c$vend[b])
            all_c$j_region[b]<-substr(all_c$nSeqCDR3[b],all_c$vend[b]+1,all_c$cd3end[b])
          }else{
            all_c$j_region[b]<-substr(all_c$nSeqCDR3[b],all_c$vend[b]+1,all_c$cd3end[b])
          }
        }
      }
    }
  }

  all_c$vregion[all_c$barcode==""] <- ""
  all_c$j_region[all_c$barcode==""] <- ""
  all_c[is.na(all_c)] <- ""
  all_c$dbegin[all_c$dbegin=="Inf"] <- ""
  all_c <- all_c[,c(1:10,12,14,13,16,15,27,11,17:26)]
  all_c$OOF <- ""
  all_c$OOF[grep('\\*|\\_', all_c$aaSeqCDR3)] <- ifelse(substr(all_c$allVHitsWithScore[grep('\\*|\\_', all_c$aaSeqCDR3)],1,4)=="TRBV","OOFB","OOFA")
  
  if(sum(all_c$OOF=="")!=length(all_c$OOF=="")){
    for(i in 1:length(unique(all_c$aaSeqCDR3[all_c$OOF!=""]))){
      all_c$OOF[all_c$aaSeqCDR3==unique(all_c$aaSeqCDR3[all_c$OOF!=""])[i]]<-paste0(all_c$OOF[all_c$aaSeqCDR3==unique(all_c$aaSeqCDR3[all_c$OOF!=""])[i]],i)}
  }
  all_c<-all_c[,c(1:5,11,28,6:10,12:16,17:27)]
  write.xlsx(all_c,"FinalOutput_MiXCR_Clonotypes.xlsx",sheetName=maj_clonotypes[a],append=T)
}


