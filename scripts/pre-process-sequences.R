######################################################################################
# format and label sequences for aligning ----
Sequences.200205<- function(dir.name= DATA)
{
  require(data.table)
  require(gdata)
  require(ape)
  
  `%notin%` <- Negate(`%in%`)
  
  #	read sequence csv data file and preprocess
  dir.name		<- '~/Box Sync/Roadmap/Data/data_220331'
  file			<- file.path(dir.name,"SHM_2201_ROADMAP_220331_tblLAB_RES.csv")
  #file.new			<- file.path(dir.name,"SHM_2201_ROADMAP_220317_tblLAB_RES.csv")
  infile.indinfo <- file.path(dir.name,'SHM_2201_ROADMAP_220331_tblBAS.csv')
  outfile			<- file.path(dir.name,"SHM_2201_ROADMAP_220331_tblLAB_seq.rda")

  verbose			<- 1	
  NA.DateRes		<- c("1911-11-11","")
  NA.Seq			<- c('')
  min.seq.len		<- 21
  proc.GeneCode	<- c(1,2)
  
  df <- read.csv(file, stringsAsFactors=FALSE)

  df <- as.data.table(df)
  setnames(df, colnames(df), toupper(colnames(df)))
  setnames(df, 'TEST_ID', 'SEQ_ID')
  #	process sequencing dates
  tmp <- df[, which(SAMPLE_D%in%NA.DateRes)]
  if(verbose) cat(paste("\nentries with missing sequence dates, n=", length(tmp)))
  #	exclude seqs with missing dates
  df <- subset(df, !SAMPLE_D%in%NA.DateRes)
  #	exclude integrase
  df <- subset(df, SEQTYPE!='IN')
  
  set(df, NULL, 'SAMPLE_D', df[,as.Date(SAMPLE_D, format="%Y-%m-%d")])
  if(verbose) cat(paste("\nrange of sequence dates is",paste(range(df[,SAMPLE_D], na.rm=1),collapse=', ')))	
 
  # remove duplicated fragments
  df <- unique(df)
  
  #	define partial pol
  df				<- df[, {
    if(length(SEQ_NUC)>2) stop('more than two sequence fragments for same ID?', SEQ_ID[1])
    z	<- ''
    if(any(SEQTYPE=='PRO'))
      z<- paste0(z,SEQ_NUC[SEQTYPE=='PRO'])
    if(any(SEQTYPE=='RT'))
      z<- paste0(z,SEQ_NUC[SEQTYPE=='RT'])				
    list(SEQ_D=SAMPLE_D[1], 						
          SEQUENCE=z)
  }, by=c('PATIENT','SEQ_ID')]

  #	define gene length
  tmp			<- df[, list(SEQ_L=nchar(SEQUENCE)), by=c('PATIENT','SEQ_ID')]	
  df			<- merge(df, tmp, by=c('PATIENT','SEQ_ID'))
  #	warning if no sequence
  if( nrow(subset(df, SEQ_L==0)) ) 
    cat('\nwarning: no sequence found for', subset(df, SEQ_L==0)[, paste(SEQ_ID, collapse=', ')])
  df			<- subset(df, SEQ_L>0)
  #	fixup seq IDs
  set(df, NULL, 'SEQ_ID', df[, gsub(' ','-',SEQ_ID)])
  
  # Distinguish Amsterdam seqs from rest of NL
  db <- read.csv(infile.indinfo, stringsAsFactors=FALSE)
  db <- as.data.table(db)
  setnames(db, colnames(db), toupper(colnames(db)))
  #ams <- db$PATIENT[db$REG_FRS_AMSTERDAM==1 | db$REG_LST_AMSTERDAM==1]
  ams <- paste0('T',db$PATIENT[db$REG_EVR_AMSTERDAM==1])
  
  #	add AMST to sequence labels
  set(df, NULL, 'SEQ_LABEL', df[, paste0('NL_',PATIENT,'_',SEQ_D,'_',SEQ_ID)])
  df$SEQ_LABEL[df$PATIENT %in% ams] <- gsub('NL_','Amst_',df$SEQ_LABEL[df$PATIENT %in% ams])
  #
  #	done reading. transform data to two objects: sequence info (data.table) and sequences (DNAbin)
  #
  ds			<- copy(df)
  set(ds, NULL, 'SEQUENCE', ds[, gsub('?','n',tolower(SEQUENCE),fixed=1)])
# Remove trail of Ns at start/end of sequences
  startswithn <- startsWith(ds$SEQUENCE,'n',ignore.case=T)
  ds$SEQUENCE[startswithn==T] <- gsub('n','-',tolower(ds$SEQUENCE[startswithn==T]),fixed=1) 
  endswithn <- endsWith(tolower(ds$SEQUENCE),'n')
  ds$SEQUENCE[endswithn==T] <- gsub('n','-',tolower(ds$SEQUENCE)[endswithn==T],fixed=1)
  save(ds, file=file.path(dir.name,"SHM_2201_ROADMAP_220331_tblLAB_seq_wSEQ.rda"))		
  
  seq			<- as.DNAbin(strsplit(ds[, SEQUENCE],''))
  names(seq)	<- ds[,SEQ_LABEL]
  write.dna(seq, file=gsub('rda','fasta',outfile), format='fasta', colsep='', nbcol=-1)
  
  ds[, SEQUENCE:=NULL]		
  #
  #	save Dutch sequences
  #	
  save(ds, file=outfile)		
  
  quit("no")	
}

# identify which sequences are new ----
find.new.seqs <- function(){
  require(ape)
  require(gdata)
  require(data.table)
  require(plyr)
  
  `%notin%` <- Negate(`%in%`)
  
  #	read sequence csv data file and preprocess
  dir.name		<- '~/Box Sync/Roadmap/Data/data_220331'
  file			<- file.path("/Users/alexb/Box Sync/Roadmap/Data/data_200821/SHM_1902_ROADMAP_200821_tblLAB_seq_wSEQ.rda")
  file.new			<- file.path(dir.name,"SHM_2201_ROADMAP_220331_tblLAB_seq_wSEQ.rda")
  infile.indinfo <- file.path(dir.name,'SHM_2201_ROADMAP_220331_tblBAS.csv')
  outfile			<- file.path(dir.name,"SHM_2201_ROADMAP_220331_tblLAB_new_seq.rda")
  
  load(file)
  df <- ds
  load(file.new)
  dn <- ds

  # count number of seqs per patient
  count <- dn[, list(N=length(SEQ_ID)),by=c('PATIENT')]
  table(count$N)
  
  df[, flag:=1]
  setnames(df,'PATIENT','PATIENT_OLD')
  # count old data
  co <- df[, list(N=length(SEQ_ID)),by=c('PATIENT_OLD')]
  table(co$N)
  
  dn <- merge(dn,subset(df,select=c('PATIENT_OLD','SEQ_D','SEQUENCE','flag')),by=c('SEQ_D','SEQUENCE'),all.x=T)

  # some patients from old data sets have new sequences but now have new PAT ID. Need to exclude those from new patients
  map <- unique(subset(dn,select=c('PATIENT','PATIENT_OLD')))
  map <- merge(map,co,by=c('PATIENT_OLD'),all.x=T)
  map <- subset(map,N<=0)
  
}

# align using virulign
#virulign '/Users/alexb/Box Sync/Roadmap/Data/K03455-pol.fasta' '/Users/alexb/Box Sync/Roadmap/Data/data_200821/SHM_2002_ROADMAP_200901_tblLAB_seq.fasta'   --exportKind GlobalAlignment   --exportAlphabet Nucleotides --exportReferenceSequence yes --nt-debug Failed > '/Users/alexb/Box Sync/Roadmap/Data/data_200821/SHM_2002_ROADMAP_200901_tblLAB_seq_aligned.fasta'
#####################################################################################
# add failed sequences to alignment ----
Add.failed.sequences.to.alignment.200205<- function(dir.name= DATA) {
require(ape)
require(gdata)
require(data.table)
require(plyr)

dir.name		<- '/Users/alexb/Box Sync/Roadmap/Data/data_200821'
file			<- file.path(dir.name,"SHM_1902_ROADMAP_200821_tblLAB_seq_aligned_trimmed.fasta")
fileall			<- file.path(dir.name,"SHM_1902_ROADMAP_200821_tblLAB_seq.fasta")
outfile			<- file.path(dir.name,"SHM_1902_ROADMAP_200821_tblLAB_seq_aligned_trimmed_complete.fasta")

seqs <- read.dna(file,format="fasta",as.matrix=T)
ids <- labels(seqs)
names(seqs) <- labels(seqs)
allseqs <- read.dna(fileall,format="fasta",as.matrix=T)
failed <- allseqs[names(allseqs) %notin% ids]
write.dna(failed, file=file.path(dir.name,"SHM_1902_ROADMAP_200821_tblLAB_seq_failed.fasta"), format='fasta', colsep='', nbcol=-1)

`%notin%` <- Negate(`%in%`)
#mafft '/Users/alexb/Box Sync/Roadmap/Data/data_200821/SHM_1902_ROADMAP_200821_tblLAB_seq_HXB2.fasta' > '/Users/alexb/Box Sync/Roadmap/Data/data_200821/SHM_1902_ROADMAP_200821_tblLAB_seq_HXB2_aligned_mafft.fasta'  

failed <- list.files("/Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223/Failed",pattern="fasta")
files <- length(failed)

# load mafft alignments
mafftalign <- read.dna("~/Box Sync/Roadmap/Data/data_200821/SHM_1902_ROADMAP_200821_tblLAB_seq_HXB2_aligned_mafft.fasta",format="fasta")

# just keep the ones which failed
seqs_mafft <- copy(mafftalign)
idsmafft <- labels(seqs_mafft)
names(seqs_mafft) <- labels(seqs_mafft)
seqs_mafft <- seqs_mafft[c(grep("HXB2",rownames(seqs_mafft)),which(rownames(seqs_mafft) %notin% rownames(seqs))),]

seqs_complete <- seq.align.based.on.common.reference(seqs,seqs_mafft, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
write.dna(seqs_complete, file=outfile, format='fasta', colsep='', nbcol=-1)
}
# Check these manually in an alignment viewer; make any manual modifications and save as alignments_DATE_Amsterdam.fasta



# ---- add failed sequences to alignment ----
Sequences.removeDRMs.200221<- function(dir.name= DATA)
{
  require(ape)
  require(big.phylo)
  dir.name		<- '~/Box Sync/Roadmap/Data/data_200821/'
  dir.name <- '/rds/general/project/ratmann_roadmap_data_analysis/live/analysis_200827/alignments'
  file <- file.path(dir.name,'SHM_1902_ROADMAP_200821_tblLAB_alignment_200908.fasta')
  outfile	<- file.path(dir.name,'SHM_1902_ROADMAP_200821_tblLAB_alignment_200908_noDRM.fasta')
  
  seq	<- read.dna(file, format='fa')
  tmp	<- which(grepl("HXB2",rownames(seq)))
  rownames(seq)[tmp]	<- 'HXB2'
  tmp	<- big.phylo:::seq.rm.drugresistance(seq)
  nodr.info	<- tmp$nodr.info
  seq	<- tmp$nodr.seq
  write.dna(seq, file= outfile, format='fasta')
  save(seq, nodr.info, file= gsub('fasta','rda',outfile))
}

# ---- remove gaps to blast ----

Sequences.no.gap.blast.200219 <- function(dir.name= DATA){
  dir.name <- '~/Box Sync/Roadmap/Data/data_200821'
  dir.name <- '/rds/general/project/ratmann_roadmap_data_analysis/live/analysis_200827/alignments'
  infile <- file.path(dir.name,'SHM_1902_ROADMAP_200821_tblLAB_alignment_200908_noDRM.fasta')
  outfile <- file.path(dir.name,'SHM_1902_ROADMAP_200821_tblLAB_alignment_200908_noDRM_nogaps.fasta')
  ds <- read.dna(infile, format='fa')
  ds			<- seq.replace(ds, code.from='-', code.to='n', verbose=0)
  write.dna(ds,outfile,format='fasta', colsep='', nbcol=-1)
}

# ---- subtype classification ----

Sequences.subtype.classification.200218 <- function(dir.name= DATA)
{
  require(ape)
  require(data.table)
  require(plyr)
  
  # Classify aligned and trimmed sequences using Comet
  # comet_osx -bs 1000 < '/Users/alexb/Box Sync/Roadmap/Data/data_200821/SHM_1902_ROADMAP_200821_tblLAB_alignment_200908_noDRM.fasta' > '/Users/alexb/Box Sync/Roadmap/Data/data_200821/SHM_1902_ROADMAP_200821_tblLAB_alignment_200908_noDRM_COMET.csv'
  
  dir.name		<- '~/Box Sync/Roadmap/Data/data_200821'
  file.c			<- file.path(dir.name,"SHM_1902_ROADMAP_200821_tblLAB_alignment_200821_noDRM_COMET.csv")
  file <- file.path(dir.name,"SHM_1902_ROADMAP_200821_tblLAB_alignment_200908_noDRM.fasta")
  outfile     <- file.path(dir.name,"sequences_to_reclassify.fasta")
  
  st		<- as.data.table(read.csv(file=file.c, sep=',', stringsAsFactors=FALSE))
  st <- st[,!'virus']
  setnames(st, c('name','subtype','bootstrap.support'), c('FASTASampleCode','SUBTYPE_C','SUBTYPE_CBS'))
  seqs <- read.dna(file,format='fasta')
  
  subtypes.c <- count(st,'SUBTYPE_C')
  subtypes.c <- subtypes.c[order(-subtypes.c$freq),]
  
  # Choose those with <35 sequences each for re-classification
  unclassified <- subtypes.c$SUBTYPE_C[subtypes.c$freq<35]
  tmp <- st[st$SUBTYPE_C %in% unclassified,]
  tokeep <- tmp$FASTASampleCode

  # Also pick out the ones which say "check" or "unassigned"
  unclassified <- c(tokeep,st$FASTASampleCode[grepl('check',st$SUBTYPE_C)])
  df <- seqs[labels(seqs) %in% unclassified,]
  write.dna(df, file=outfile, format='fasta', colsep='', nbcol=-1)
  
  unclassified <- c(st$FASTASampleCode[grepl('unassigned',st$SUBTYPE_C)])
  df <- seqs[labels(seqs) %in% unclassified,]
  write.dna(df, file=outfile, format='fasta', colsep='', nbcol=-1)
  
  # Upload file to REGA subtyping tool and obtain fasta output
  infile <- file.path(dir.name,'REGA_classifications_3.csv')
  outfile.c <- file.path(dir.name,'REGA_classifications.rds')
  outfile <- file.path(dir.name,'Complete_subtype_classifications.rds')
  r <- read.delim(infile,sep=",")
  saveRDS(r,file=outfile.c)
  subtypes.r <- count(r,'assignment')
  subtypes.r <- subtypes.r[order(-subtypes.r$freq),]
  r.a <- r[,c(1,3)]

  # Update subtype classification with REGA results
  setnames(r.a, 'name', 'FASTASampleCode')
  setnames(r.a, 'assignment', 'SUBTYPE_R')
  st <- merge(st,r.a,by='FASTASampleCode',all.x=T)
  # Replace original subtypes with those from REGA
  # Make classifications consistent e.g. HIV-1 Subtype B = B (12_BF vs. CRF 12_BF)
  
  st$SUBTYPE <- st$SUBTYPE_C
  st$SUBTYPE[st$SUBTYPE_R!='NA'] <- 'Other'
  
  # update with Rega subtype for those with <100% BS support
  set(st, st[, which(SUBTYPE_R=='HIV-1 Subtype B' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', 'B')
  set(st, st[, which(SUBTYPE_R=='HIV-1 Subtype C' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', 'C')
  set(st, st[, which(SUBTYPE_R=='HIV-1 Subtype D' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', 'D')
  set(st, st[, which(SUBTYPE_R=='HIV-1 Subtype F (F1)' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', 'F1')
  set(st, st[, which(SUBTYPE_R=='HIV-1 Subtype G' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', 'G')
  set(st, st[, which(SUBTYPE_R=='HIV-1 Subtype H' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', 'H')
  set(st, st[, which(SUBTYPE_R=='HIV-1 Subtype K' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', 'K')
  set(st, st[, which(SUBTYPE_R=='HIV-1 Subtype A (A1)' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', 'A1')
  set(st, st[, which(SUBTYPE_R=='HIV-1 Subtype A (01_AE)' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', '01_AE')
  set(st, st[, which(SUBTYPE_R=='HIV-1 CRF 02_AG' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', '02_AG')
  set(st, st[, which(SUBTYPE_R=='Recombinant of A1, G' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', 'A1-G')
  set(st, st[, which(SUBTYPE_R=='Recombinant of G, A1' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', 'G-A1')
  set(st, st[, which(SUBTYPE_R=='Recombinant of B, F1' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', 'B-F1')
  set(st, st[, which(SUBTYPE_R=='HIV-1 CRF 06_CPX' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', '06_cpx')
  set(st, st[, which(SUBTYPE_R=='Recombinant of B, A1' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', 'B-A1')
  set(st, st[, which(SUBTYPE_R=='HIV-1 CRF 07_BC' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', '07_BC')
  set(st, st[, which(SUBTYPE_R=='HIV-1 Subtype A (01_AE)' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', '01_AE')
  set(st, st[, which(SUBTYPE_R=='HIV-1 CRF 01_AE' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', '01_AE')
  set(st, st[, which(SUBTYPE_R=='Recombinant of 02_AG, A1' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', '02_AG-A1')
  set(st, st[, which(SUBTYPE_R=='HIV-1 CRF 12_BF' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', '12_BF')
  set(st, st[, which(SUBTYPE_R=='Recombinant of B, C' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', 'B-C')
  set(st, st[, which(SUBTYPE_R=='HIV-1 Subtype B-like' & (SUBTYPE_CBS!=100 | is.na(SUBTYPE_CBS)))], 'SUBTYPE', 'B')
  
  set(st, st[, which(SUBTYPE_C=='A1 (check for 02_AG)' & SUBTYPE_R=='HIV-1 CRF 02_AG')], 'SUBTYPE', '02_AG')
  set(st, st[, which(SUBTYPE_C=='A1 (check for 02_AG)' & SUBTYPE_R=='HIV-1 CRF 02_AG-like A1')], 'SUBTYPE', '02_AG')
  set(st, st[, which(SUBTYPE_C=='A1 (check for 02_AG)' & SUBTYPE_R=='HIV-1 CRF 02_AG-like')], 'SUBTYPE', '02_AG')
  set(st, st[, which(SUBTYPE_C=='A1 (check for 02_AG)' & SUBTYPE_R=='Recombinant of A1, 02_AG')], 'SUBTYPE', '02_AG-A1')
  set(st, st[, which(SUBTYPE_C=='A1 (check for 02_AG)' & SUBTYPE_R=='HIV-1 Subtype A (A1), potential recombinant')], 'SUBTYPE', 'A1')
  set(st, st[, which(SUBTYPE_C=='A1 (check for 02_AG)' & SUBTYPE_R=='HIV-1 Subtype A (A1)')], 'SUBTYPE', 'A1')
  set(st, st[, which(SUBTYPE_C=='A1 (check for 02_AG)' & SUBTYPE_R=='HIV-1 Subtype A (A1)-like')], 'SUBTYPE', 'A1')
  
  set(st, st[, which(SUBTYPE_C=='G (check for 02_AG)' & SUBTYPE_R=='HIV-1 CRF 02_AG')], 'SUBTYPE', '02_AG')
  set(st, st[, which(SUBTYPE_C=='G (check for 02_AG)' & SUBTYPE_R=='HIV-1 Subtype G')], 'SUBTYPE', 'G')
  set(st, st[, which(SUBTYPE_C=='G (check for 02_AG)' & SUBTYPE_R=='HIV-1 Subtype G-like')], 'SUBTYPE', 'G')
  set(st, st[, which(SUBTYPE_C=='G (check for 11_cpx)' & SUBTYPE_R=='HIV-1 Subtype G-like')], 'SUBTYPE', 'G')
  set(st, st[, which(SUBTYPE_C=='G (check for 02_AG)' & SUBTYPE_R=='HIV-1 Subtype G (02_AG)')], 'SUBTYPE', '02_AG')
  set(st, st[, which(SUBTYPE_C=='G (check for 02_AG)' & SUBTYPE_R=='HIV-1 CRF 02_AG-like')], 'SUBTYPE', '02_AG')
  set(st, st[, which(SUBTYPE_C=='G (check for 02_AG)' & SUBTYPE_R=='HIV-1 CRF 06_CPX-like')], 'SUBTYPE', '06_cpx')
  
  set(st, st[, which(SUBTYPE_C=='F1' & SUBTYPE_R=='HIV-1 Subtype F (F1)')], 'SUBTYPE', 'F1')
  set(st, st[, which(SUBTYPE_C=='F1 (check for 71_BF1)' & SUBTYPE_R=='HIV-1 Subtype F (F1)-like')], 'SUBTYPE', 'F1')
  set(st, st[, which(SUBTYPE_C=='72_BF1' & SUBTYPE_R=='HIV-1 Subtype F (F1)-like')], 'SUBTYPE', 'F1')
  
  set(st, st[, which(SUBTYPE_C=='B (check for 29_BF)' & SUBTYPE_R=='HIV-1 Subtype B-like')], 'SUBTYPE', 'B')
  set(st, st[, which(SUBTYPE_C=='B (check for 90_BF1)' & SUBTYPE_R=='HIV-1 Subtype B-like')], 'SUBTYPE', 'B')
  set(st, st[, which(SUBTYPE_C=='B (check for 29_BF)' & SUBTYPE_R=='HIV-1 Subtype B-like')], 'SUBTYPE', 'B')
  set(st, st[, which(SUBTYPE_C=='56_cpx' & SUBTYPE_R=='HIV-1 Subtype B-like')], 'SUBTYPE', 'B')
  set(st, st[, which(SUBTYPE_C=='44_BF' & SUBTYPE_R=='HIV-1 Subtype B-like')], 'SUBTYPE', 'B')
  set(st, st[, which(SUBTYPE_C=='38_BF1' & SUBTYPE_R=='HIV-1 Subtype B-like')], 'SUBTYPE', 'B')
  
  set(st, st[, which(SUBTYPE_C=='C (check for 92_C2U)' & SUBTYPE_R=='HIV-1 Subtype C-like')], 'SUBTYPE', 'C')
  
  set(st, st[, which(SUBTYPE_C=='01_AE (check for 15_01B)' & SUBTYPE_R=='HIV-1 CRF 01_AE')], 'SUBTYPE', '01_AE')
  
  set(st, st[, which(SUBTYPE_C=='unassigned_1, 01_AE-B' & SUBTYPE_R=='HIV-1 CRF 01_AE')], 'SUBTYPE', '01_AE')
  set(st, st[, which(SUBTYPE_C=='unassigned_1' & SUBTYPE_R=='HIV-1 Subtype A (A1)-like')], 'SUBTYPE', 'A1')
  set(st, st[, which(SUBTYPE_C=='unassigned_1' & SUBTYPE_R=='HIV-1 Subtype G-like')], 'SUBTYPE', 'G')
  set(st, st[, which(SUBTYPE_C=='unassigned_1' & SUBTYPE_R=='HIV-1 Subtype C-like')], 'SUBTYPE', 'C')
  set(st, st[, which(SUBTYPE_C=='unassigned_1' & SUBTYPE_R=='HIV-1 Subtype B-like')], 'SUBTYPE', 'B')
  set(st, st[, which(SUBTYPE_C=='unassigned_2' & SUBTYPE_R=='HIV-1 CRF 02_AG')], 'SUBTYPE', '02_AG')
  set(st, st[, which(SUBTYPE_C=='unassigned_2' & SUBTYPE_R=='HIV-1 Subtype G-like')], 'SUBTYPE', 'G')
  set(st, st[, which(SUBTYPE_C=='unassigned_2' & SUBTYPE_R=='HIV-1 CRF 06_CPX-like')], 'SUBTYPE', '06_cpx')
  
  set(st, st[, which(SUBTYPE_C=='HIV-1 Subtype A (A1), potential recombinant' & SUBTYPE_R=='A1 (check for 02_AG)')], 'SUBTYPE', 'Other')
  set(st, st[, which(SUBTYPE_C=='Recombinant of F1, A1' & SUBTYPE_R=='A1 (check for 02_AG)')], 'SUBTYPE', 'Other')
  set(st, st[, which(SUBTYPE_C=='HIV-1 Subtype A (A1)-like' & SUBTYPE_R=='A1 (check for 02_AG)')], 'SUBTYPE', 'Other')
  set(st, st[, which(SUBTYPE_R=='HIV-1 CRF 09_cpx')], 'SUBTYPE_C', 'Other')
  set(st, st[, which(SUBTYPE_R=='Check the report')], 'SUBTYPE_C', 'Other')
  set(st, st[, which(grepl('unassigned',SUBTYPE_C) & SUBTYPE_R=='Check the report')], 'SUBTYPE', 'Other')

  set(st, st[, which(SUBTYPE=='02_AG-A1')], 'SUBTYPE', 'Other')
  set(st, st[, which(SUBTYPE=='G-A1')], 'SUBTYPE', 'Other')
  set(st, st[, which(SUBTYPE=='B-F1')], 'SUBTYPE', 'Other')
  set(st, st[, which(SUBTYPE=='B-A1')], 'SUBTYPE', 'Other')
  set(st, st[, which(SUBTYPE=='07_BC')], 'SUBTYPE', 'Other')
  set(st, st[, which(SUBTYPE=='K')], 'SUBTYPE', 'Other')
  set(st, st[, which(SUBTYPE=='12_BF')], 'SUBTYPE', 'Other')
  set(st, st[, which(SUBTYPE=='B-C')], 'SUBTYPE', 'Other')
  set(st, st[, which(SUBTYPE=='H')], 'SUBTYPE', 'Other')
  
  # check counts of subtypes
  subtypes.cr <- plyr::count(st,'SUBTYPE')
  subtypes.cr <- subtypes.cr[order(-subtypes.cr$freq),]
  subtypes.cr
  
  # check most common for Amsterdam seqs only
  ams <- subset(st,grepl("Amst",st$FASTASampleCode))
  subtypes.cr <- plyr::count(ams,'SUBTYPE')
  subtypes.cr <- subtypes.cr[order(-subtypes.cr$freq),]
  subtypes.cr
  
  saveRDS(st,file=outfile)
}

# ---- process LANL seqs ----
Sequences.LANL.preprocess.200217 <- function(dir.name= DATA)
{
  require(ape)
  require(gdata)
  require(data.table)
  
  dir.name		<- '/Users/alexb/Documents/Roadmap/Data/blast_db'
  infile <- file.path(dir.name,'data_200205_LANL_raw.fasta')
  outfile <- file.path(dir.name,'data_200205_LANL_maxlen.fasta')

  ds <- read.dna(infile, format='fa')
  ds <- as.character(ds)
  ds.len <- max(sapply(ds, length))
  
  # Add gaps to end of sequences to make same length
  seqs <- length(ds) 
  for(i in 1:seqs)
  {
    if(i %% 1e3 == 0)
      print(paste0('done: ', i))
    z <- ds[[i]]
    z <- paste0( c(z, rep('-', ds.len - length(z) ) ), collapse='' ) 
    cat('>',names(ds)[i],'\n',z,'\n', file=outfile, append=TRUE)	
  }
  
  infile			  <- file.path(dir.name,"data_200205_LANL_maxlen.fasta")
  outfileNL			<- file.path(dir.name,"data_200205_LANL_noNL.fasta")
  
  ds <- read.dna(infile,format="fasta")
  
  # Remove dutch sequences
  ids <- labels(ds)
  dutch <- grepl('NL.',ids,fixed=1)
  df <- ds[!dutch,]
  df <- as.list(df)
  
  # Takes a while to save
  write.dna(df, file=outfileNL, format='fasta', colsep='', nbcol=-1)
  
  dir.name		<- '/Users/alexb/Documents/Roadmap/Data/blast_db'
  infile			<- file.path(dir.name,"data_200205_LANL_noNL.fasta")
  outfile.lanl <- file.path(dir.name,'data_2002_LANL_alignment.fasta')
  outfile.nogaps <- file.path(dir.name,'data_2002_LANL_alignment_nogaps.fasta')
  outfile.blast <- file.path(dir.name,'data_2002_LANL_alignment_blast.fasta')
  outfile.nodrm.blast <- file.path(dir.name,'data_2002_LANL_alignment_noDRM_blast.fasta')
  
  # Strip the gaps to make the file smaller
  s			  <- read.dna(infile, format='fa')

  nbatch <- ceiling(ncol(s)/1e3)
  for(j in 1:nbatch)
  {
    if (j==1) {s2      <- seq.strip.gap(s[,seq.int(9001, 9394)],strip.pc=0.99998, gap.chars='-')
    } else {s2      <- seq.strip.gap(s[,seq.int((nbatch-j)*1e3+1, (nbatch-j)*1e3 + 1e3)],strip.pc=0.9998, gap.chars='-')}
    batchfile <- file.path('/Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223/batchfiles',paste(j,".fasta",sep=""))
    write.dna(s2,file=batchfile,format='fasta', colsep='', nbcol=-1)
  }
  # concatenate the fasta files
  j <- 10
  batchfile <- file.path('/Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223/batchfiles',paste(j,".fasta",sep=""))
  z <- read.dna(batchfile,format='fasta') 
  for(j in 9:1)
  {
    batchfile <- file.path('/Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223/batchfiles',paste(j,".fasta",sep=""))
    s2 <- read.dna(batchfile,format='fasta') 
    z <- cbind(z,s2)
  }
  write.dna(tmp,file=outfile.nogaps, format='fasta', colsep='', nbcol=-1) 
  
  # Just keep columns in HXB2
  z <- read.dna(outfile.nogaps,format='fa')
  ids <- labels(z)
  hxb2 <- grepl('K03455',ids,fixed=1)
  z	<- as.character(z)
  nogap <- which(z[hxb2,]!='-')
  z <- as.DNAbin(z[,nogap,drop=FALSE])
  
  z <- read.dna(outfile.lanl,format='fa')
  # Remove any gaps from end of taxa labels
  rownames(dbs) <- gsub('([ ]+)$','',rownames(dbs))
  # Remove taxa which did not align well
  toremove <- c(
  "01_AE.CN.2009.8059.MK286989",
  "02_AG.KW.2016.65991_16.KX155646",
  "02H.AO.1993.93AOHDC251.KU310618",
  "07_BC.CN.2009.967.KJ820361",
  "07_BC.CN.2009.967.KP235147",
  "B.US.2009.86GW09P5C11.MK385498",
  "C.ZA.2014.Pt4.1_S1426.MK643780",
  "BG.ES.2003.X138-6.EU074779",
  "-.BW.-.074-A-110-1-4_w72_264.MK458028",
  "A1.RW.2007.R880FPL_d157_A5.KP223798",
  "19_cpx.CU.2011.CUPL1101.JX194180",
  "B.US.2011.72GW11P4E5.MK385378",
  "B.US.2014.37GW14P3C12.MK384862",
  "B.US.2012.2302-PBCM-37.KY778461",
  "C.ET.2008.ET147.KU319541",
  "C.ET.2008.ET155.KU319546",
  "A1O.FR.2006.BCF212.KY359380",
  "O.FR.1992.VAU.AF407418",
  "-.FR.2013.RBF235.MN337382",
  "13_cpx.GL.2004.GRL100.AM285274",
  "BF1.IT.-.621A.AY855645",
  "BF1.AR.-.13.AY140111",
  "HU.PT.2003.PTCBR_49.HM135484",
  "02_AG.FR.-.201.HM035625",
  "B.TH.2010.R75018.KC962530",
  "C.ET.2009.ETH-G-5648.KF026147",
  "18_cpx.CU.2009.CU1163-09.JQ585277",
  "B.CU.2010.CUPL1035.JN000007",
  "-.BW.-.074-C-310-1-2_w65_448.MK458212",
  "O.SN.1998.SE42HALD.AJ300450",
  "O.US.-.DEOXXUS001.KF859744",
  "O.DE.-.DEOXXDE004.KF859742",
  "O.ES.2001.Read25_HIV_GroupO.KX228804",
  "O.ES.-.ESP1.U97171",
  "C.ZA.2007.ZA_2007_DGM_3.JN176288",
  "N.CM.1999.YBF116.AJ564923",
  "04_cpx.GR.2004.BP00053_SUP_LH01.JN687664",
  "O.SN.1998.SE42HALD.AJ300450",
  "B.CU.2014.14CU006.KR914677",
  "DO.FR.2008.RBF208.GQ351296",
  "O.ES.-.DEOXXES001.KF859743",
  "O.DE.-.DEOXXDE004.KF859742",
  "O.ES.2001.Read25_HIV_GroupO.KX228804",
  "O.BE.1987.ANT70.L20587"
  )
  tokeep <- z[!(labels(z) %in% toremove),]
  write.dna(tokeep,file=outfile.lanl, format='fasta', colsep='', nbcol=-1) 
  
  # Remove DRMs
  lfile <- file.path('/Users/alexb/Documents/Roadmap/Data/blast_db/data_200305_LANL_alignment.fasta')
  loutfile <- file.path('/Users/alexb/Documents/Roadmap/Data/blast_db/data_200305_LANL_alignment_noDRM.fasta')
  seq					<- read.dna(lfile, format='fa')
  tmp					<- which(grepl("K03455",rownames(seq)))
  rownames(seq)[tmp]	<- 'HXB2'
  tmp					<- seq.rm.drugresistance(seq)
  nodr.info			<- tmp$nodr.info
  seq					<- tmp$nodr.seq
  write.dna(seq, file= loutfile, format='fasta')
  save(seq, nodr.info, file= gsub('fasta','rda',loutfile))
  
  # keep first sequence per patient
  labs <- data.table(TAXA=rownames(seq))
  labs$PID <- gsub('([A-Za-z0-9_]+).([A-Za-z]+).([0-9]+).([A-Za-z0-9_]+).([A-Za-z0-9]+)','\\4',labs$TAXA_L)
  db$PID2 <- gsub('([A-Za-z0-9_]+).([A-Za-z]+).([0-9]+).([A-Za-z0-9_]+).([A-Za-z0-9]+)','\\4',db$TAXA_L)
  
  
  labs$ACC <- gsub('([A-Za-z0-9_-]+).([A-Za-z-]+).([0-9-]+).([A-Za-z0-9_-]+).([A-Za-z0-9-]+)','\\5',labs$TAXA_L)
  labs$YEAR <- gsub('([A-Za-z0-9_-]+).([A-Za-z-]+).([0-9-]+).([A-Za-z0-9_-]+).([A-Za-z0-9-]+)','\\3',labs$TAXA_L)
  
  setkey(labs, PID, YEAR)
  first	<- labs[, list(SEQUENCE=TAXA[1]), by='PID']
  aseq <- aseq[labels(aseq) %in% c(hxb2lab,first$SEQUENCE),]
  
  # Replace gaps with 'n' for BLAST db
  seq			<- seq.replace(seq, code.from='-', code.to='n', verbose=0)
  write.dna(tokeep,file=outfile.blast, format='fasta', colsep='', nbcol=-1) 
  write.dna(tokeep,file=outfile.nodrm.blast, format='fasta', colsep='', nbcol=-1) 
}

# ---- make tree of LANL ref seqs to find closest subtype/CRF for outgroup ----

Subgroup.outgroups.similar.subtypes<- function(dir.name= DATA)
{
  require(ape)
  dir.name		<- '/Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223'
  rfile <- file.path(dir.name,'HIV1_SUBTYPE_REF_2010_2253-3870_DNA.fasta')
  ds <- read.dna(rfile,format='fa')
  # Draw tree for ref sequences to find similar subtypes
  tree <- nj(dist)
  class(tree)
  # Just keep part of taxa label with subtype
  #tree$tip.label <- substr(tree$tip.label,1,regexpr("\\.",tree$tip.label)-1)
  tree$tip.label <- substr(tree$tip.label,1,10)
  plot(tree, cex = 0.4)
}

# ---- add closest LANL sequences and make subtype alignments ----

Sequences.addcloseLANL.200218<- function(dir.name= DATA)
{
  require(plyr)
  require(dplyr)
  require(tidyr)
  require(ape)
  require(data.table)
  
  dir.name		<- '~/Box Sync/Roadmap/Data/data_200821'
  file.c			<- file.path(dir.name,"SHM_1902_ROADMAP_200821_tblLAB_alignment_200908_noDRM_COMET.csv")
  file <- file.path(dir.name,"SHM_1902_ROADMAP_200821_tblLAB_alignment_200908_noDRM.fasta")
  
  # Create blast DB
  #makeblastdb -in /Users/alexb/Documents/Roadmap/Data/blast_db/data_200305_LANL.fasta -dbtype nucl -title LANL_polany_n35000_strippedn -out /Users/alexb/Documents/Roadmap/Data/blast_d/LANL_polany_n35000_strippedn.db
  
  # Run blast query
  #blastn -query /Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223/SHM_1902_ROADMAP_200821_tblLAB_alignment_200821_noDRM_nogaps.fasta -db /Users/alexb/Documents/Roadmap/Data/blast_db/LANL_polany_n82000_strippedn.db -out /Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223/SHM_1902_ROADMAP_200821_tblLAB_alignment_200821_LANL_BLAST_20.txt -max_target_seqs 20 -outfmt 6
  
  dir.name		<- '~/Box Sync/Roadmap/Data/data_200821'
  bfile <- file.path(dir.name,'SHM_1902_ROADMAP_200821_tblLAB_alignment_200821_LANL_BLAST_20.txt')
  lfile <- file.path(dir.name,'data_200305_LANL_alignment_noDRM.fasta')
  afile <- file.path(dir.name,'SHM_1902_ROADMAP_200821_tblLAB_alignment_200821_noDRM.fasta')
  sfile <- file.path(dir.name,'Complete_subtype_classifications.rds')
  infile.seqinfo			<- file.path(dir.name,"SHM_1902_ROADMAP_200821_tblLAB_seq.rda")
  infile.indinfo <- file.path(dir.name,'SHM_1902_ROADMAP_191223_tblBAS.csv')
  outfile.matches      <- file.path(dir.name,'ROADMAP_200821_All_Taxa_withsubtype.rda')
  outfile.seqs <- file.path(dir.name,'ROADMAP_200821_All_Sequences.fasta')
  outfile.trimseqs <- file.path(dir.name,'ROADMAP_200821_All_Sequences_LANL_aligned.fasta')

  st <- readRDS(sfile)
  dl <- read.dna(lfile,format='fa')
  db		<- as.data.table(read.delim(bfile, sep='\t',stringsAsFactors=FALSE, header=FALSE))		
  setnames(db, paste('V',1:12,sep=''),c('FASTASampleCode','TAXA_L','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'))
  # Sort (ascending)
  setkey(db, FASTASampleCode, pident)
  db$rowid <- rownames(db)
  
  # Obtain 20 closest LANL sequences for each Ams sequence
  tmp <- db[, {
      z<- seq(length(pident),max(1,length(pident)-19),-1) 
      list(TAXA_L=TAXA_L[z], pident=pident[z], rowid=rowid[z])
    }, by='FASTASampleCode']
  # Merge back
  db <- merge(db, tmp, by=c('FASTASampleCode','TAXA_L','pident','rowid'))
  db <- subset(db, select = -c(rowid))
  
  # Get subtypes for each Amsterdam sequence  
  db	<- merge(db, subset(st, select=c(FASTASampleCode, SUBTYPE)), by='FASTASampleCode', all.x=T)
  
  # Get subtype of each LANL match
  db <- as.data.table(db)
  db[, SUBTYPE_L := substr(TAXA_L,1,regexpr("\\.",TAXA_L)-1)]

  # Re-order by closest match
  setkey(db, FASTASampleCode, pident)
  
  # Get subtype of closest match
  best	<- db[, list(subtype.lc=SUBTYPE_L[length(SUBTYPE_L)]), by='FASTASampleCode']
  ds		<- merge(db, best, by=c('FASTASampleCode')) 
  save(ds,file=outfile.matches)  
  
  # Keep the unique ones only
  dbc			<- unique(subset(ds, select=TAXA_L))
  dbs			<- read.dna(lfile, format='fa')
  rownames(dbs) <- gsub('([ ]+)$','',rownames(dbs))
  hxb2lab <- dbc[,TAXA_L][grepl('K03455',dbc[,TAXA_L],fixed=1)]
  # Rename HXB2 with Taxa label for matching
  rownames(dbs)[rownames(dbs)=='HXB2']	<- hxb2lab
  dba <- dbs[labels(dbs) %in% c(dbc[, TAXA_L],hxb2lab),]
  
  # Make sure HXB2 is in the file
  ids <- rownames(dba)
  hxb2 <- grepl('K03455',ids,fixed=1)
  hxb2lab <- ids[hxb2]
  
  # Add Dutch sequences to LANL ones
  aseq <- read.dna(afile,format='fa')
  ids <- rownames(aseq)
  hxb2 <- grepl('HXB2',ids,fixed=1)
  rownames(aseq)[hxb2] <- hxb2lab

  # Only keep first sequence per patient
  load(infile.seqinfo)
  dseq <- ds
  dseq <- as.data.table(dseq)
  setkey(dseq, PATIENT, SEQ_D)
  first	<- dseq[, list(SEQUENCE=SEQ_LABEL[1]), by='PATIENT']
  aseq <- aseq[labels(aseq) %in% c(hxb2lab,first$SEQUENCE),]
  
  # Add LANL seqs
 seqs.all <- seq.align.based.on.common.reference(aseq, dba, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')

  write.dna(seqs.all,file=outfile.seqs, format='fasta', colsep='', nbcol=-1) 
  seqs.all <- seqs.all[,1:1302]
  write.dna(seqs.all,file=outfile.trimseqs, format='fasta', colsep='', nbcol=-1) 
  
  # Read back in
  sfile <- file.path(dir.name,'ROADMAP_200821_All_Sequences_LANL_aligned.fasta')
  seqs.all <- read.dna(outfile.trimseqs,format='fa')
  outfile.matches <- file.path(dir.name,'ROADMAP_200821_All_Taxa_withsubtype.rda')
  load(outfile.matches)
    
  # Amend subtypes and remove singletons (after checking trees)
  # Remove singletons and very short sequences
  toremove <- c(
    
  )
  ds <- ds[!toremove,]
  
  # Read in outgroup reference sequences
  rfile <- file.path(dir.name,'HIV1_SUBTYPE_REF_2010_2253-3870_DNA.fasta')
  rs <- read.dna(rfile,format='fa')
  rownames(rs) <- gsub('Ref','Outgroup',rownames(rs))
  ids <- labels(rs)
  
  #	select closest sequences for subtype B
  outfile.b <- file.path(dir.name,'alignments','ROADMAP_200407_Sequences_LANL_aligned_noDRM_subtype_B_wOutgroup.fasta')
  db.l			<- unique(subset(ds, SUBTYPE_L=='B', TAXA_L))
  db.a			<- unique(subset(ds, SUBTYPE=='B', FASTASampleCode))
  db.b <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode]),]
  st.b <- ids[c(grep('K03455',ids,fixed=1),grep('.28_BF.',ids,fixed=1))] 
  sub.b <- rs[labels(rs) %in% st.b,]
  scn				<- seq.align.based.on.common.reference(db.b, sub.b, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
  scn <- scn[,1:1302]
  write.dna(scn, file=outfile.b, format='fa', colsep='', nbcol=-1)
  
  # Subtype 02_AG
  outfile.02ag <- file.path(dir.name,'alignments','ROADMAP_200407_Sequences_LANL_aligned_noDRM_subtype_02AG_wOutgroup.fasta')
  db.l			<- unique(subset(ds, SUBTYPE_L=='02_AG', TAXA_L))
  db.a			<- unique(subset(ds, SUBTYPE=='02_AG', FASTASampleCode))
  db.02ag <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],hxb2lab),]
  st.02ag <- ids[c(grep('K03455',ids,fixed=1),grep('.36_cpx.',ids,fixed=1))] 
  sub.02ag <- rs[labels(rs) %in% st.02ag,]
  scn				<- seq.align.based.on.common.reference(db.02ag, sub.02ag, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
  scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
  scn <- scn[,1:1302]
  write.dna(scn, file=outfile.02ag, format='fa', colsep='', nbcol=-1)
  
  # Subtype C
  outfile.c <- file.path(dir.name,'alignments','ROADMAP_200407_Sequences_LANL_aligned_noDRM_subtype_C_wOutgroup.fasta')
  db.l			<- unique(subset(ds, SUBTYPE_L=='C', TAXA_L))
  db.a			<- unique(subset(ds, SUBTYPE=='C', FASTASampleCode))
  db.c <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],hxb2lab),]
  st.c <- ids[c(grep('K03455',ids,fixed=1),grep('.31_BC.',ids,fixed=1))] 
  sub.c <- rs[labels(rs) %in% st.c,]
  scn				<- seq.align.based.on.common.reference(db.c, sub.c, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
  scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
  scn <- scn[,1:1302]
  write.dna(scn, file=outfile.c, format='fa', colsep='', nbcol=-1)
  
  # Subtype A1
  outfile.a1 <- file.path(dir.name,'alignments','ROADMAP_200407_Sequences_LANL_aligned_noDRM_subtype_A1_wOutgroup.fasta')
  db.l			<- unique(subset(ds, SUBTYPE_L=='A1', TAXA_L))
  db.a			<- unique(subset(ds, SUBTYPE=='A1', FASTASampleCode))
  db.a1 <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],hxb2lab),]
  st.a1 <- ids[c(grep('K03455',ids,fixed=1),grep('.35_AD.',ids,fixed=1))] 
  sub.a1 <- rs[labels(rs) %in% st.a1,]
  scn				<- seq.align.based.on.common.reference(db.a1, sub.a1, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
  scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
  scn <- scn[,1:1302]
  write.dna(scn, file=outfile.a1, format='fa', colsep='', nbcol=-1)
  
  # Subtype 01_AE
  outfile.01ae <- file.path(dir.name,'alignments','ROADMAP_200407_Sequences_LANL_aligned_noDRM_subtype_01AE_wOutgroup.fasta')
  db.l			<- unique(subset(ds, SUBTYPE_L=='01_AE', TAXA_L))
  db.a			<- unique(subset(ds, SUBTYPE=='01_AE', FASTASampleCode))
  db.01ae <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],hxb2lab),]
  st.01ae <- ids[c(grep('K03455',ids,fixed=1),grep('.15_01B.',ids,fixed=1))] 
  sub.01ae <- rs[labels(rs) %in% st.a1,]
  scn				<- seq.align.based.on.common.reference(db.01ae, sub.01ae, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
  scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
  scn <- scn[,1:1302]
  write.dna(scn, file=outfile.01ae, format='fa', colsep='', nbcol=-1)
  
  # Subtype G
  outfile.g <- file.path(dir.name,'alignments','ROADMAP_200407_Sequences_LANL_aligned_noDRM_subtype_G_wOutgroup.fasta')
  db.l			<- unique(subset(ds, SUBTYPE_L=='G', TAXA_L))
  db.a			<- unique(subset(ds, SUBTYPE=='G', FASTASampleCode))
  db.g <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],hxb2lab),]
  st.g <- ids[c(grep('K03455',ids,fixed=1),grep('.14_BG.',ids,fixed=1))] 
  sub.g <- rs[labels(rs) %in% st.g,]
  scn				<- seq.align.based.on.common.reference(db.g, sub.g, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
  scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
  scn <- scn[,1:1302]
  write.dna(scn, file=outfile.g, format='fa', colsep='', nbcol=-1)
  
  # Subtype D
  outfile.d <- file.path(dir.name,'alignments','ROADMAP_200407_Sequences_LANL_aligned_noDRM_subtype_D_wOutgroup.fasta')
  db.l			<- unique(subset(ds, SUBTYPE_L=='D', TAXA_L))
  db.a			<- unique(subset(ds, SUBTYPE=='D' , FASTASampleCode))
  db.d <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],hxb2lab),]
  st.d <- ids[c(grep('K03455',ids,fixed=1),grep('.19_cpx.',ids,fixed=1))] 
  sub.d <- rs[labels(rs) %in% st.d,]
  scn				<- seq.align.based.on.common.reference(db.d, sub.d, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
  scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
  scn <- scn[,1:1302]
  write.dna(scn, file=outfile.d, format='fa', colsep='', nbcol=-1)
  
  # Subtype 06_cpx
  outfile.06cpx <- file.path(dir.name,'alignments','ROADMAP_200407_Sequences_LANL_aligned_noDRM_subtype_06cpx_wOutgroup.fasta')
  db.l			<- unique(subset(ds, SUBTYPE_L=='06_cpx', TAXA_L))
  db.a			<- unique(subset(ds, SUBTYPE=='06_cpx', FASTASampleCode))
  db.06cpx <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],hxb2lab),]
  st.06cpx <- ids[c(grep('K03455',ids,fixed=1),grep('.45_cpx.',ids,fixed=1))]
  sub.06cpx <- rs[labels(rs) %in% st.06cpx,]
  scn				<- seq.align.based.on.common.reference(db.06cpx, sub.06cpx, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
  scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
  scn <- scn[,1:1302]
  write.dna(scn, file=outfile.06cpx, format='fa', colsep='', nbcol=-1)
  
}
