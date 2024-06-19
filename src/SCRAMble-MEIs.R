#!/usr/bin/env Rscript
# analyzes clusters(reads in cluster files) for MEIs
# parses command line options
library(stringr)
library(Biostrings)
library(IRanges)
suppressPackageStartupMessages(library(dplyr))
suppressMessages(require(Rsamtools))
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, c("-o", "--out-file"), type="character", dest="out.file", help="full path to output file")
parser <- add_option(parser, c("-c", "--cluster-file"), type="character", dest="cluster.file",  help="full path to cluster file")
parser <- add_option(parser, c("-n", "--nCluster"), type="integer", dest="nCluster", default=5, help="min cluster size to analyze [default %default]")
parser <- add_option(parser, c("-m", "--mei-score"),  type="integer", default=50, dest="mei.score", help="min mei alignment score to call [default %default]")
parser <- add_option(parser, c("--poly-a-frac"),  type="numeric", default=0.75, dest="polyAFrac", help="fraction of clipped length for calling polyA tail in MEIs [default %default]")
parser <- add_option(parser, c("--pct-align"),  type="integer", default=90, dest="pctAlign", help="percent alignment of clipped read for calling deletionss [default %default]")
parser <- add_option(parser, c("--poly-a-dist"),  type="integer", default=100, dest="polyAdist", help="how far from MEI to look for polyA tail [default %default]")
parser <- add_option(parser, c("--mei-refs"), type="character", default="MEI_consensus_seqs.fa", dest="mei.ref.file", help="full path to MEI reference file (fasta format) [default %default]")
parser <- add_option(parser, c("-r", "--ref"), type="character", default=NULL, dest="ref", help="reference file (fasta format) [default %default]")
opt = parse_args(parser)
## init globals from arguments
outFile         = opt[["out.file"]]
clusterFile     = opt[["cluster.file"]]
nCluster        = opt[["nCluster"]]
meiScore        = opt[["mei.score"]]
mei.refs        = opt[["mei.ref.file"]]
polyAFrac       = opt[["polyAFrac"]]
pctAlign        = opt[["pctAlign"]]
polyAdist       = opt[["polyAdist"]]
blastRef        = opt[["ref"]]
# print arguments
cat("-----------------------------------analyzes clusters(reads in cluster files) for MEIs-----------------------------------\n")
cat("Running sample:", clusterFile, "\n")
objects = ls()
objects = objects[-match(c("opt", "parser"),objects)]
cat("Running MEIs with options:\n")
for(i in 1:length(objects)){
  cat(objects[i],":", get(objects[i]), "\n")
}
#
# Given a string or vector of strings, returns the reverse compelement
#' rc('ACGGGACTACTTAACG')
rc = function(string){
  # fxn to take reverse complement of a single string or vector of strings
  r.string = tolower(reverse(string))
  rc.string = gsub("a", "T", r.string)
  rc.string = gsub("t", "A", rc.string)
  rc.string = gsub("c", "G", rc.string)
  rc.string = gsub("g", "C", rc.string)
  return(tolower(rc.string))
}
# convert fasta reference MEI sequences into mobilome R object
make.mobilome = function(mei.refs){
  refs = read.delim(mei.refs, as.is=T, header=F)
  mobilome = list()
  indices = which(grepl("^>", refs[,1]))
  for(i in 1:length(indices)){
    seq = paste0(refs[(indices[i]+1):ifelse(i<length(indices),(indices[i+1]-1),nrow(refs)),1], collapse="")
    mobilome = c(mobilome, list(tolower(seq)))
    names(mobilome)[i] = gsub(">", "", as.character(refs[indices[i],1]))
  }
  return(mobilome)
}
# Takes the clipped read consensus (forward sequence and its reverse complement) and aligns to the mobilome.
# Alignment quality statistics are returned. Filtering thresholds are used to determine likely MEIs.
aligner = function(my.fa, ref=mobilome, rc=F){
  results = {}
  for(i in 1:length(ref)){
    alignment=pairwiseAlignment(pattern = tolower(my.fa), subject = ref[[i]], type="local")
    alignment.df = data.frame(alignment_len=nchar(alignment),
                              alignment_score=score(alignment),
                              percent_identity = pid(alignment), stringsAsFactors = F)
    alignment.df$percent_clipped_read_aligned = (alignment.df$alignment_len * 100) / nchar(my.fa)
    alignment.df$starts = alignment@subject@range@start
    alignment.df$stops = alignment.df$starts + alignment@subject@range@width - 1
    alignment.df$MEI_Family = names(ref)[i]
    if(i==1){
      results = alignment.df
    }else{
      results = rbind.data.frame(results, alignment.df)
    }
  }
  results$rc = rc
  return(results)
}
# Make pretty output table
pick.winners = function(df, pct_read_aligned=75, meiScore=50){
  cols = c("coord","MEI_Family", "alignment_score", "percent_clipped_read_aligned", "percent_identity", "starts", "stops", "rc")
  my.winners = df[df$alignment_score >= meiScore & df$percent_clipped_read_aligned >= pct_read_aligned,]
  if(nrow(my.winners)==0){
    winners = data.frame(Insertion=character(),
                         MEI_Family=character(),
                         Insertion_Direction=character(),
                         ClippedReads_In_Cluster=character(),
                         Alignment_Score=character(),
                         Alignment_Percent_Length=character(),
                         Alignment_Percent_Identity=character(),
                         ClippedSequence=character(),
                         Clipped_Side=character(),
                         Start_In_MEI=character(),
                         Stop_In_MEI=character(),
                         RNAME= character(),
                         clipped_pos = character(),
                         stringsAsFactors = F)
  }else{
    winners = data.frame(Insertion=my.winners$coord,
                         MEI_Family=my.winners$MEI_Family,
                         Insertion_Direction=ifelse(!my.winners$rc, "Plus", "Minus"),
                         ClippedReads_In_Cluster=my.winners$counts,
                         Alignment_Score=my.winners$alignment_score,
                         Alignment_Percent_Length=my.winners$percent_clipped_read_aligned,
                         Alignment_Percent_Identity=my.winners$percent_identity,
                         ClippedSequence=toupper(my.winners$clipped.consensus),
                         Clipped_Side=my.winners$clipped,
                         Start_In_MEI=my.winners$starts,
                         Stop_In_MEI=my.winners$stops,
                         RNAME= my.winners$RNAME,
                         clipped_pos = my.winners$clipped_pos,
                         stringsAsFactors = F)
  }
  return(winners)
}

do.meis = function(all,  refs, polyAFrac=0.5, meiScore=50,
                   pctAlign=70, polyAdist=200){
  mobilome = make.mobilome(refs)
  df.all = all[all$counts >= nCluster , ] # INCLUDE SIMPLE SEQUENCE FOR MEI SEARCH
  df.all$clipped.consensus.rc = rc(df.all$clipped.consensus)
  alignments_fwd = aligner(my.fa=df.all$clipped.consensus, ref=mobilome)
  alignments_rev = aligner(my.fa=df.all$clipped.consensus.rc, ref=mobilome, rc=T)
  df.aligned = rbind.data.frame(data.frame(df.all, alignments_fwd,stringsAsFactors = F, row.names=NULL),
                                data.frame(df.all, alignments_rev, stringsAsFactors = F, row.names=NULL) )
  # Make pretty output table
  df.aligned$coord = paste(df.all$RNAME, df.all$clipped_pos, sep=":")
  winners = as.data.frame(pick.winners(df.aligned,  pct_read_aligned=pctAlign, meiScore=meiScore) )
  cat("Sample had", nrow(winners), "MEI(s)\n")
  # IDENTIFY POLY-A TAIL AND TSD
  if(nrow(winners)==0){
    winners = data.frame(winners, polyA_Position=character(), polyA_Seq=character(),
                         polyA_SupportingReads=character(), TSD=character(), TSD_length=character())
  }else{
    winners$polyA_Position = "None Found"
    winners$polyA_Seq = "None Found"
    winners$polyA_SupportingReads = "None Found"
    winners$TSD = "None Found"
    winners$TSD_length = "NA"
  }
  if(nrow(winners)>0){
    for(i in 1:nrow(winners)){
      if(is.na(winners$Insertion_Direction[i]) | winners$Insertion_Direction[i] == "Plus"){
        # then look for left clipped with poly-A upstream
        polyA = all[all$clipped=="left" & all$RNAME == winners$RNAME[i] & all$clipped_pos >= winners$clipped_pos[i]-polyAdist & all$clipped_pos <= winners$clipped_pos[i],]
        if(nrow(polyA)==0){ next }
        polyA$A_count = str_count(polyA$clipped.consensus, "a")
        polyA$A_frac = polyA$A_count/nchar(polyA$clipped.consensus)
        if( any(!is.na(polyA$A_frac) & polyA$A_frac >= polyAFrac)){
          # grab the one closest to the MEI
          polyA=polyA[order(polyA$clipped_pos, decreasing=T),]
          polyA = polyA[polyA$A_frac>=polyAFrac,][1,]
          winners$TSD[i] = toupper(substr(polyA$anchored.consensus, start=1, stop=(winners$clipped_pos[i] - polyA$clipped_pos ) ))
          winners$TSD_length[i] = winners$clipped_pos[i] - polyA$clipped_pos
          if((winners$clipped_pos[i] - polyA$clipped_pos ) > nchar(polyA$anchored.consensus)){
            winners$TSD[i] = paste0(winners$TSD[i], "<TRUNCATED>",  collapse="")
          }
        }else{next()}
      }else{
        # then look for right clipped with poly-T downstream
        polyA = all[all$clipped=="right" & all$RNAME == winners$RNAME[i] & all$clipped_pos <= winners$clipped_pos[i]+polyAdist & all$clipped_pos >= winners$clipped_pos[i],]
        if(nrow(polyA)==0){ next }
        polyA$A_count = str_count(polyA$clipped.consensus, "t")
        polyA$A_frac = polyA$A_count/nchar(polyA$clipped.consensus)
        if(any(!is.na(polyA$A_frac) & polyA$A_frac >= polyAFrac)){
          # grab the one closest to the MEI
          polyA=polyA[order(polyA$clipped_pos, decreasing=F),]
          polyA = polyA[polyA$A_frac>=polyAFrac,][1,]
          winners$TSD_length[i] = polyA$clipped_pos - winners$clipped_pos[i]
          if(1 > (nchar(polyA$anchored.consensus) - (polyA$clipped_pos - winners$clipped_pos[i]))){
            winners$TSD[i] = paste0("<TRUNCATED>",toupper(polyA$anchored.consensus), collapse="")
          }else{
            winners$TSD[i] = toupper(substr(polyA$anchored.consensus, start= (nchar(polyA$anchored.consensus) - (polyA$clipped_pos - winners$clipped_pos[i]) + 1), stop=nchar(polyA$anchored.consensus)) )
          }
        } else{next()}
      }
      winners$polyA_Position[i] = polyA$clipped_pos
      winners$polyA_Seq[i] = toupper(polyA$clipped.consensus)
      winners$polyA_SupportingReads[i] = polyA$counts
    }
  }
  winners = winners[,-(match(c("RNAME", "clipped_pos"), names(winners)) )]
  winners = winners[order(winners$Insertion, winners$Alignment_Score, decreasing=T),]
  winners = winners[!duplicated(winners$Insertion),]
  return(winners)
}
# convert to VCF
mei.write = function(SampleID, winners, fa, blastRef=None, outFile=None){
  if (is.null(winners) | nrow(winners) == 0 | missing(fa)) return(NULL)
  header = c('##fileformat=VCFv4.2',
             paste('##reference=', blastRef, sep=''),
             paste('##contig=<ID=', names(fa), '>', sep=""),
             '##FILTER=<ID=PASS,Description="All filters passed">',
             '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
             '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
             '##INFO=<ID=END,Number=.,Type=Integer,Description="End position for structural variants">',
             '##INFO=<ID=Strands,Number=.,Type=String,Description="Strand orientation of the adjacency">',
             '##INFO=<ID=ClippedSequence,Number=.,Type=String,Description="ClippedSequence">',
             '##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">',
             '##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">',
             '##ALT=<ID=INS:ME:SVA,Description="Insertion of SVA element">',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             '##FORMAT=<ID=CR,Number=1,Type=Integer,Description="Number of clipped reads supporting the variant">'
             )
  fixed = data.frame('#CHROM' =  gsub("(.*):(\\d*)$", "\\1", winners$Insertion),
                    'POS' = as.integer(gsub("(.*):(\\d*)$", "\\2", winners$Insertion)),
                    'ID' = paste(sub(':', '-', winners$Insertion), toupper(winners$MEI_Family), winners$Insertion_Direction, sep="-"),
                    'FILTER' = 'PASS',
                    'ALT' = paste('<INS:ME:', toupper(winners$MEI_Family), '>', sep=''),
                    'QUAL' = winners$Alignment_Score,
                    'FORMAT' = 'GT:CR',
                    'SampleID' = paste('./.', winners$ClippedReads_In_Cluster, sep=':'),
                    'Strands' = ifelse(winners$Insertion_Direction == 'Plus', "+", "-"),
                    'ClippedSequence' = winners$ClippedSequence,
                    stringsAsFactors = F,check.names = F)
  fixed$INFO = paste(paste('Strands', fixed$Strands, sep='='), paste('ClippedSequence', fixed$ClippedSequence, sep='='), sep=';')
  fixed$REF = sapply(1:nrow(fixed), function(i) as.vector(subseq(fa[fixed[i, '#CHROM']], start=fixed$POS[i], end=fixed$POS[i])))
  if(nrow(fixed) > 0){
    fixed = fixed[, c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SampleID')]
  }else{
    fixed = data.frame('#CHROM' = character(),
               'POS' = character(),
               'ID' = character(),
               'REF' = character(),
               'ALT' = character(),
               'QUAL' = character(),
               'FILTER' = character(),
               'INFO' = character(),
               'FORMAT' = character(),
               'SampleID' = character(),
               check.names = F)
  }
  colnames(fixed) <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', SampleID)
  desired_order <- c(names(fa))
  # desired_order <- c('chrM', paste0('chr',as.character(1:22)), 'chrX', 'chrY')
  fixed <- fixed %>% filter(`#CHROM` %in% desired_order) %>% arrange(match(`#CHROM`, desired_order), POS)
  # suppressWarnings(write.table(winners, sub('.vcf', '.txt', outFile), row.names=F, col.names=T, quote=F, append=F, sep='\t'))
  suppressWarnings(write.table(header, outFile, row.names=F, col.names=F, quote=F, append=F, sep='\t'))
  suppressWarnings(write.table(fixed, outFile, row.names=F, col.names=T, quote=F, append=T, sep='\t'))
}

if(is.null(blastRef)){
  warning("A reference .fa file is required for writing results to VCF")
}else{
  # Load reference fasta
  fa = getSeq(open(FaFile(blastRef)))
  names(fa) = gsub("\\s.*", "", names(fa))
  #
  all = read.delim(clusterFile, as.is=T, header=F, col.names=c("coord", "clipped", "counts", "clipped.consensus", "anchored.consensus"))
  chrom = gsub("(.*):\\d*$", "\\1", all$coord)
  pos = pmin(pmax(as.integer(gsub(".*:", "", all$coord)), 1), width(fa[chrom]))
  all$coord = paste(chrom, pos, sep=":")
  all$clipped_pos = as.integer(gsub(".*:", "", all$coord))
  all$RNAME = gsub("(.*):\\d*$", "\\1", all$coord)
  all$rname_clippedPos_Orientation_ReadSide = paste(all$RNAME, all$clipped_pos, all$clipped, sep="_")
  SampleID = c(sub("\\..*", "", basename(clusterFile)))
  mei.winners = cbind(SampleID=SampleID, do.meis(all=all, meiScore=meiScore, refs=mei.refs, polyAFrac=polyAFrac, pctAlign=pctAlign, polyAdist=polyAdist))
  mei.write(SampleID, mei.winners, fa, blastRef=blastRef, outFile=outFile)
  cat("Done analyzing MEIs, writing VCF file to ", outFile, "\n")
}