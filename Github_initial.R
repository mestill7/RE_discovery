# module load gnu-4.7.2/bedtools-2.22.0
# module load R-3.4.1
R
library(GenomicFeatures)
library(GenomicRanges)
library(GenomicAlignments)
library(Rsamtools)
library(biovizBase)
library(rtracklayer)
# source("http://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Hsapiens.UCSC.hg38")
require("BSgenome.Hsapiens.UCSC.hg38")
require(BSgenome)
library("bedr")
library("biomaRt")

setwd(mycurrentdirectory)
sample.names.file="./SRR_Acc_List.txt"
samples.location="./align"
out.directory="./RE_discovery_v28"
gtf.file= "~/hg38/gencode.v28.primary.annotation.gtf"
mybsgenome="BSgenome.Hsapiens.UCSC.hg38"

##Helper functions
##Used for processing bam files into bigwigs and calculating number of aligned reads.
getbw <- function(i,samples.i=samples,sampleloc.i=sampleloc,suffix.i=suffix,needchrom.i=needchrom,od.i=od){
	# i=1
	# return(i)
	# return(missing(samples))
	x=samples.i[i]
	# return(x)
	message(paste("Processing sample",x))
	bamfile <- paste0(sampleloc.i,"/",x,suffix.i)
	# return(bamfile)
	alignment <- readGAlignments(bamfile)
	# return(alignment)
	reads_coverage <- coverage(alignment)
	if(!(is.null(needchrom.i))){
			if(length(needchrom.i)==length(needchrom.i[needchrom.i%in%names(reads_coverage)])){
				message("All requested chromosomes are in your bamfile")
				alignment <- alignment[seqnames(alignment)%in%needchrom.i] ##cut down alignments to needchrom
				reads_coverage <- reads_coverage[needchrom.i]
				export.bw(reads_coverage, con = paste0(od.i,"/",x,"_BigWig.bw"))
			}
			if(length(names(reads_coverage)[names(reads_coverage)%in%needchrom.i])<length(needchrom.i)){
				message("Your bamfile is missing one or more of the requested chromosomes")
				alignment <- alignment[seqnames(alignment)%in%needchrom.i]
				reads_coverage <- reads_coverage[names(reads_coverage)%in%needchrom.i]
				export.bw(reads_coverage, con = paste0(od.i,"/",x,"_BigWig.bw"))
			}
	}else{
		message("You are exporting all chromosomes of your bamfile")
		export.bw(reads_coverage, con = paste0(od.i,"/",x,"_BigWig.bw"))
	}
	return(c(x,length(alignment)))
}

##Used for retrieving coverage from bigwigs
getCVG <- function(i,od.i=od,samples.i=samples,sp.i=sp,exons.i=exons){ ##Coverage function with genome checks
	# i=1
ar=import.bw(paste(od.i,"/",samples.i,"_BigWig.bw",sep="")[i],as = c("RleList"))
message("Comparing bigwig to input genome")
if(any(!(seqlevels(sp.i)%in%seqlevels(ar)))){
	warning("Genome chromosomes not in bigwig")
	if(length(seqlevels(ar))>length(seqlevels(sp.i))){
		message("Reducing your bigwig to input genome")
		ar <- ar[seqlevels(ar)%in%seqlevels(sp.i)]
	}
	if(length(seqlevels(ar))<length(seqlevels(sp.i))){
	message("Less chromosomes in bigwig than input genome. Reducing input genome...")
	sp.i <- sp.i[seqnames(sp.i) %in% seqlevels(ar)]
	seqlevels(sp.i) <- seqlevels(sp.i)[seqlevels(sp.i)%in%seqlevels(ar)]
	}
}
ar <- ar[seqlevels(sp.i)] ##ensure ar and sp.i are in same order
res <- binnedAverage(sp.i, ar, "mean_cvg") ##Coverage across bins
res <- res[-(queryHits(findOverlaps(res, exons.i,ignore.strand=TRUE))),] ##Remove tiles that overlap exons
return(res)
}

##Used for concatenating ensembl gene IDs
concatENS <- function(x){
    myinter <- inter[inter$ID==x,]
    myinter.uniq <- myinter[!(duplicated(myinter$ID))]
    myinter.uniq$gene_id <- paste(unique(myinter$grl_name),collapse=",")
    return(myinter.uniq)
}

##Used for extending exons into near-exon REs
extendCVG <- function(i){
# i=1
newblips <- list()
newexons <- list()
ar <- import.bw(paste0(od,"/",samples,"_BigWig.bw")[i],as = c("RleList"))
message(paste("Processing sample",i,"of",length(samples)))
message("Reducing bigwig to match chromosomes with exons of interest")
inboth <- intersect(seqlevels(nn.subj),names(ar))
seqlevels(nn.subj, pruning.mode="coarse") <- inboth
ar <- ar[inboth]
seqlevels(nn.qr, pruning.mode="coarse") <- inboth
res <- binnedAverage(nn.subj, ar, "mean_cvg") ##Coverage across bins
message("Processing exons and adjacent elements")
for(x in 1:length(res)){
	# x=1
	myex <- res[x,] ##exonic region
	myqr <- nn.qr[x,] ##novel region
	myqr.bin <- unlist(tile(myqr,width=10), recursive = TRUE, use.names = TRUE) ##cut novel region into bins
	myit <- myqr.bin@seqnames@lengths ##how many novel bins do we have?
if(myit==1){ ## Only one 10 bp blip 
	nn.s <- distanceToNearest(myex,myqr.bin,ignore.strand=TRUE)
	myex2 <- myex[queryHits(nn.s)] ##create the exon to append
	myqr2 <- myqr.bin[subjectHits(nn.s)] ##create
	res2 <- binnedAverage(myqr2, ar, "mean_cvg")
	exon_cvg <- myex@elementMetadata@listData$mean_cvg
	novel_cvg <- res2@elementMetadata@listData$mean_cvg
	ratio <- abs(exon_cvg-novel_cvg)/exon_cvg
	if(exon_cvg==0){ ##exon isn't expressed, so no need to extend it.
		myint <- myex2
		myqr.bin <- myqr.bin
	}else if(ratio<0.5){
		myint <- reduce(c(myex2,myqr2,ignore.mcols=TRUE),min.gapwidth=50,ignore.strand=TRUE)
		myqr.bin <- setdiff(myqr.bin,myint,ignore.strand=TRUE) 
	}else{
		myint <- myex2
		myqr.bin <- myqr.bin
	}
}

if(myit>1){ ##multiple 10 bp blips
	nn.s <- distanceToNearest(myex,myqr.bin,ignore.strand=TRUE)
	myex2 <- myex[queryHits(nn.s)] ##create the exon to append
	myqr2 <- myqr.bin[subjectHits(nn.s)] ##create
	res2 <- binnedAverage(myqr2, ar, "mean_cvg")
	exon_cvg <- myex@elementMetadata@listData$mean_cvg
	novel_cvg <- res2@elementMetadata@listData$mean_cvg
	ratio <- abs(exon_cvg-novel_cvg)/exon_cvg
	# myit <- myqr.bin@seqnames@lengths ##how many novel bins do we have?
	if(exon_cvg==0){ ##exon isn't expressed, so no need to extend it.
		myint <- myex2
		myqr.bin <- myqr.bin
	}else if(ratio>=0.5){
		myint <- myex2
		myqr.bin <- myqr.bin
		# print("sup")
	}else{ ##if ratio <0.05
		repeat{
			myint <- reduce(c(myex2,myqr2,ignore.mcols=TRUE),min.gapwidth=50,ignore.strand=TRUE) ##combines the exon & novel regions
			# new_cvg <- binnedAverage(myint, ar, "mean_cvg")@elementMetadata@listData$mean_cvg
			myqr.bin <- setdiff(myqr.bin,myint,ignore.strand=TRUE)##Remove the added blip from the myqr.bin
			myqr.bin <- unlist(tile(myqr.bin,width=10), recursive = TRUE, use.names = TRUE)##tile the blip
			nn.s <- distanceToNearest(myint,myqr.bin,ignore.strand=TRUE)
			myqr2 <- myqr.bin[subjectHits(nn.s)] ##Grab the blip of interest
			res2 <- binnedAverage(myqr2, ar, "mean_cvg") ##Coverage of blip
			exon_cvg <- binnedAverage(myint, ar, "mean_cvg")@elementMetadata@listData$mean_cvg
			novel_cvg <- res2@elementMetadata@listData$mean_cvg
			myex2 <- myint ##refresh the exonic region for next loop
			ratio <- abs(exon_cvg-novel_cvg)/exon_cvg
			myit <- myqr.bin@seqnames@lengths ##how many novel bins do we have?
		if(ratio>=0.5){
			break
		}##end break statement
		if(ratio<0.5 & myit==1){
			myint <- reduce(c(myex2,myqr2,ignore.mcols=TRUE),min.gapwidth=50,ignore.strand=TRUE) ##combines the exon & novel regions
			# new_cvg <- binnedAverage(myint, ar, "mean_cvg")@elementMetadata@listData$mean_cvg
			myqr.bin <- setdiff(myqr.bin,myint,ignore.strand=TRUE)##Remove the added blip from the myqr.bin
			break
		}##end one tile structure
		}##end repeat structure
	}##end if loop
}##end multiple 10 bp blips loop
newexons[[x]] <- myint
newblips[[x]] <- myqr.bin

}##Complete element loop
#get the info of what exon the new exons came from
newexons2 <- do.call("c",lapply(newexons,function(x) c(x,ignore.mcols=TRUE)))#drop the metadata from each individual grange and combine to one granges object
newexon.anno <- as.data.frame(res,row.names = 1:length(res))[,6:7]
newexons2 <- as.data.frame(newexons2,row.names=1:length(newexons))
newexons2$ID <- paste(newexons2$seqnames,newexons2$start,newexons2$end,sep="_")
combexon <- cbind(newexons2,newexon.anno)

##combine the blips
newblips2 <- do.call("c",lapply(newblips,function(x) c(x,ignore.mcols=TRUE)))#drop the metadata from each individual grange and combine to one granges object
newblips2 <- reduce(newblips2) ##don't need nearby exon info yet, still need to combine the results from different samples
return(list(combexon,newblips2))
} ##complete coverage function

##Reformatting function 

getex <- function(i){
	# i=1
	inter <- do.call(rbind,lapply(exlist,function(x) return(x[i,]))) #i is row number
	inter2 <- reduce(as(inter,"GRanges"),ignore.strand=TRUE)
	inter2 <- cbind(as.data.frame(inter2),inter[1,7:8])
	inter2 <- inter2[,c(1:3,6:7,5)]
	return(inter2)
}
##End helper functions


generatebw <- function(needchrom,sample.names=sample.names.file,sampleloc=samples.location,od=out.directory,suffix="_subsample.bam"){
library(GenomicAlignments)
library(rtracklayer)
samples <- read.table(sample.names,as.is=T,sep="\t")[,1]
environment(getbw) <- environment()
res <- data.frame(t(sapply(1:length(samples),getbw)),stringsAsFactors=TRUE)
# return(res)
res$X2 <- as.numeric(as.character(res$X2))
colnames(res) <- c("sample","read_count")
write.table(res,paste0(od,"/alignment_sizes.txt"),col.names=TRUE,row.names=FALSE,sep="\t")
}
generatebw(c(paste("chr",1:22,sep=""),"chrX", "chrY", "chrM"))

binGenome <- function(mygenome,genversion="hg38",species="Hsapiens"){
	library(mygenome,character.only=TRUE)
	message("Currently reading chosen genome")
	genome <- get(species)
	gen <- Seqinfo(seqnames=c(seqnames(genome)),
	             seqlengths=c(seqlengths(genome)),
	             isCircular=c(isCircular(genome)),
	             genome=genversion)
	message("Currently binning the genome")
	sp <- tileGenome(gen, tilewidth=10, cut.last.tile.in.chrom=FALSE) ##split the genome into 10 bp tiles
	sp <- unlist(sp, recursive = TRUE, use.names = TRUE)
	return(sp)
}
sp.whole <- binGenome(mybsgenome)

deriveExons <- function(needchrom,sp=sp.whole,gtf=gtf.file){
	message("Currently reading gene annotation")
	txdb <- makeTxDbFromGFF(gtf, format="gtf") ##To make TxDb from external GTF
	exons <- flatGrl(exonsBy(txdb,by=c("tx"),use.names=TRUE)) ##Turn exons into Granges object
	if(any(seqlevels(sp)!=seqlevels(exons))){
		warning("Your genome has different chromosome names than your annotation, was that on purpose?")
		tmpchr <- seqlevels(exons)[seqlevels(exons)%in%seqlevels(sp)]
		message("Removing annotation chromsomes that are NOT in your genome")
		exons <- exons[seqnames(exons) %in% tmpchr]
		seqlevels(exons) <- tmpchr
		if(!(is.null(needchrom))){
			if(any(!(needchrom%in%seqnames(sp)))){
				message("One or more requested chromsomes are not in your genome. Removing...") ##commands cutting down to input chromosomes
				needchrom <- needchrom[needchrom%in%seqnames(sp)]
				sp <- sp[seqnames(sp) %in% needchrom]
				seqlevels(sp) <- needchrom 
				exons <- exons[seqnames(exons) %in% needchrom]
				seqlevels(exons) <- seqlevels(exons)[seqlevels(exons)%in%needchrom]
			}
		}
		if(is.null(needchrom)){
			message("You have not pre-selected any chromosomes of interest. Modify \"needchrom\" object if desired.")
		}
		message(paste("You currently have",length(seqlevels(sp)),"chromosomes in your genome. If more or less than desired, modify your command."))
	}
	return(list(sp,exons))
}
exon.out <- deriveExons(c(paste("chr",1:22,sep=""),"chrX", "chrY", "chrM"))

# ar.test <- import.bw("/wsu/home/eh/eh97/eh9735/mott-ngs/IBS/dbp/test/RE_discovery/5060v2_seqH5_BigWig.bw")

findRE <- function(needchrom,sample.names=sample.names.file,sampleloc=samples.location,od=out.directory,gtf=gtf.file,mu_size=2.5,iteration_size=2,sp=exon.out[[1]],exons=exon.out[[2]]){
samples <- read.table(sample.names,as.is=T,sep="\t")[,1]
# return(samples)
ar.set <- read.table(paste0(od,"/alignment_sizes.txt"),as.is=T,header=TRUE,sep="\t")
if(all(ar.set$sample==samples)){
	mu <- ar.set$read_count*(mu_size*1/10^6)
}else{
	stop("Your read counts are out of order.")
}
environment(getCVG) <- environment()
myit <- seq(1,length(samples),by=iteration_size)
	for (itvar in myit){
		# itvar=19
		if(itvar!=max(seq(1,length(samples),by=iteration_size))){
			iteration=c(itvar:(itvar+(iteration_size-1)))
			print(itvar)
			# iteration=c(241:246)
			res_list <- lapply(iteration,getCVG)
			conCOV <- function(x){return(x$mean_cvg)}
			res2 <- do.call(cbind,(lapply(res_list,conCOV)))
			colnames(res2) <- samples[iteration]
		    res2 <- data.frame(res2,check.names=FALSE)
			good <- do.call(cbind,lapply(samples[iteration],function(x) res2[,x] >= mu[x])) ##which bins are gte the library size normalized parameter?
			rownames(good) <- rownames(res2)
			good2 <- good[apply(good,1,function(x) length(x[x=="TRUE"])>=(length(x)/2)),]
			good3 <- cbind(as.data.frame(res_list[[1]])[rownames(good2),1:5],good2) ##combine good2 with genomic region
			good3 <- as(good3,"GRanges")
			good4 <- reduce(good3,min.gapwidth=100L)
			write.table(as.data.frame(good4)[,1:3],paste(od,"/blips_iteration_",itvar,".bed",sep=""))
			}else{
				iteration=c(itvar:length(samples))
				print(itvar)
				res_list <- lapply(iteration,getCVG)
				conCOV <- function(x){return(x$mean_cvg)}
				res2 <- do.call(cbind,(lapply(res_list,conCOV)))
				colnames(res2) <- samples[iteration]
			    res2 <- data.frame(res2,check.names=FALSE)
				good <- do.call(cbind,lapply(samples[iteration],function(x) res2[,x] >= mu[x])) ##which bins are gte the library size normalized parameter?
				rownames(good) <- rownames(res2)
				good2 <- good[apply(good,1,function(x) length(x[x=="TRUE"])>=(length(x)/2)),]
				if(length(iteration)>1){
					good3 <- cbind(as.data.frame(res_list[[1]])[rownames(good2),1:5],good2) ##combine good2 with genomic region
					good3 <- as(good3,"GRanges")
					good4 <- reduce(good3,min.gapwidth=100L)
					write.table(as.data.frame(good4)[,1:3],paste(od,"/blips_iteration_",itvar,".bed",sep=""))
					}else{
						good3 <- cbind(as.data.frame(res_list[[1]])[names(good2),1:5],good2) ##combine good2 with genomic region
						good3 <- as(good3,"GRanges")
						good4 <- reduce(good3,min.gapwidth=100L)
						write.table(as.data.frame(good4)[,1:3],paste(od,"/blips_iteration_",itvar,".bed",sep=""))	
					}
			}
	}
}
# findRE(c("chr8","chrY", "chrM"))
findRE(c(paste("chr",1:22,sep=""),"chrX", "chrY", "chrM"))


combineRE <- function(od,sample.names=sample.names.file){
	files <- paste0(od,"/blips_iteration_",seq(1,length(samples),by=iteration_size),".bed")
	readfs <- function(x){return(read.table(x))}
	bps <- lapply(files,readfs)
	mer <- do.call(rbind,bps)
	mer <- as(mer,"GRanges")
	mer <- reduce(mer, min.gapwidth=50L)
	df <- as.data.frame(mer)
	df$misc <- "BLIP"
	df$name <- paste("BLIP_",rownames(df),sep="")
	df <- df[,c("seqnames","start","end","name","misc","strand")]
	write.table(df,paste0(od,"/Blips_merged.bed"),sep="\t",quote=F,col.names =F,row.names = F)
}
combineRE("/wsu/home/eh/eh97/eh9735/mott-ngs/IBS/serono/RE_discovery")

prepareExons <- function(needchrom,gtf=gtf.file,sample.names=sample.names.file,od=out.directory){
samples <- read.table(sample.names,as.is=T,sep="\t")[,1]
message("Currently reading gene annotation")
txdb <- makeTxDbFromGFF(gtf, format="gtf") ##To make TxDb from external GTF
txdb_exons <- exonsBy(txdb,by=c("gene")) ##Get the exons locations
res <- as.data.frame(flatGrl(txdb_exons))##Turn exons into Granges object
res <- res[,c("seqnames","start","end","grl_name","exon_name","strand")] ##order dataframe into bed format
##Clean up whitespaces
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
res$seqnames <- sapply(res$seqnames,trim)
res$grl_name <- sapply(res$grl_name,trim)
res$exon_name <- sapply(res$exon_name,trim)
res$strand <- sapply(res$strand,trim)
colnames(res) <- c("chr","start","end","grl_name","exon_name","strand")
##This code is designed to work with both non-"chr" and "chr" chromosomes
if(is.null(needchrom)){
    a.valid <- is.valid.region(res,check.chr = FALSE)
    res.valid <- res[a.valid,]
    is.sorted <- is.sorted.region(res.valid,check.chr = FALSE) ##FALSE
    # res.sort.natural <- bedr.sort.region(res.valid, method = "natural",check.chr = FALSE)
    res.sort.natural <- bedr.sort.region(res.valid, method = "lexicographical",check.chr = FALSE,check.merge=FALSE)
    is.merged <- is.merged.region(res.sort.natural,check.chr = FALSE)
    res.merge.stranded <- bedr.merge.region(res.sort.natural,stratify.by="strand",check.chr = FALSE)
    ##also need to remove duplicates on plus and minus strand
    res.merge.unique <- res.merge.stranded[!(duplicated(res.merge.stranded[,1:3])),]  ##removes 20 exons
}
if(!(is.null(needchrom))){ ##if needchrom is present
    res <- res[res$chr%in%needchrom,]
    a.valid <- is.valid.region(res,check.chr = FALSE)
    res.valid <- res[a.valid,]
    is.sorted <- is.sorted.region(res.valid,check.chr = FALSE) ##FALSE
    # res.sort.natural <- bedr.sort.region(res.valid, method = "natural",check.chr = FALSE)
    res.sort.natural <- bedr.sort.region(res.valid, method = "lexicographical",check.chr = FALSE,check.merge=FALSE)
    is.merged <- is.merged.region(res.sort.natural,check.chr = FALSE)
    res.merge.stranded <- bedr.merge.region(res.sort.natural,stratify.by="strand",check.chr = FALSE)
    ##also need to remove duplicates on plus and minus strand
    res.merge.unique <- res.merge.stranded[!(duplicated(res.merge.stranded[,1:3])),]  ##removes 20 exons
}
##need to add on the "exon" and strand columns
res.merge.unique$class <- "EXON"
res.merge.unique$strand <- unlist(lapply(strsplit(rownames(res.merge.unique),".", fixed = TRUE),function(x) return(x[1])))
write.table(res.merge.unique,paste0(od,"/gtf_exons.bed"),quote=F,col.names =F,row.names = F,sep="\t") 
}

prepareExons(c(paste("chr",1:22,sep=""),"chrX", "chrY", "chrM"))



##As this section uses biomaRt, you must provide the biomart and dataset needed for the "useMart" function
##several commands use ensembl IDs. manually modify the "bm" variable if necessary for your organism.
##Get locations of REs: are they intergenic or intronic
annotateRE <- function(gtf,od=out.directory,mymart="ensembl",mydataset="hsapiens_gene_ensembl"){
ensembl=useMart(biomart=mymart, dataset=mydataset)
message("Reading in annotation file, used for finding introns")
txdb <- makeTxDbFromGFF(gtf, format="gtf") ##To make TxDb from external GTF
message("Processing and renaming introns")

introns <- flatGrl(intronsByTranscript(txdb,use.names=TRUE))
##modify introns to change ENST to ENSG in this section
introns$grl_short <- unlist(lapply(strsplit(as.character(introns$grl_name),"\\."),function(x) return(x[1])))
bm <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","external_gene_name"),filters='ensembl_transcript_id', values=as.character(unique(introns$grl_short)),mart=ensembl)
if(is.null(bm) | nrow(bm)<5){
	stop("biomaRt retrieval seems to have failed.")
}
introns$grl_name <- bm[match(introns$grl_short,bm[,2]),1]
introns$grl_short <- NULL
message("Reading in pre-collapsed exons.")
exon <- read.table(paste0(od,"/gtf_exons.bed"),as.is=T,sep="\t")
colnames(exon) <- c("seqnames","start","end","grl_name","class","strand")
exon <- as(exon,"GRanges")
message("Reading in RNA elements.")
blip <- read.table(paste0(od,"/Blips_merged_hg38.bed"),sep="\t")
colnames(blip) <- c("seqnames","start","end","grl_name","class","strand")
blip <- as(blip,"GRanges")
introns.b <- setdiff(introns, exon, ignore.strand=TRUE) ##remove any exonic areas from the introns
introns.c <- introns.b[queryHits(findOverlaps(introns.b, introns,ignore.strand=TRUE)),] ##Retrieve all overlaps of *exon removed* introns with regular introns
introns.named <- introns[subjectHits(findOverlaps(introns.b, introns,ignore.strand=TRUE)),] ##Retrieve overlaps of regular introns with *exon removed* introns, to get gene names
inter <- introns.c
inter$grl_name <- introns.named$grl_name
inter$ID <- paste(seqnames(inter),start(inter),end(inter),sep="_")
##To save computation time, we preselect for just the introns that overlap a blip
##Computation with all human introns (~277k unique locations) takes ~ 5 hours
## ~ 1 minute per 1000 unique locations
hits <- findOverlaps(inter,blip,ignore.strand=TRUE)
inter <- inter[queryHits(hits)] 
uniqID <- unique(inter$ID)
message(paste("Collapsing introns. This step will take approximately",length(uniqID)/1000,"minutes."))
environment(concatENS) <- environment()
system.time(res <- lapply(uniqID,concatENS)) ##This takes ~5 hours for ~277k introns
message("Patience is a virtue. Introns successfully collapsed.")
introns.d <- do.call("c",res)

##Now overlap the blips with the introns 
message("Overlapping novel elements with introns")
hits= findOverlaps(blip, introns.d,type="within",ignore.strand=TRUE)
int.bl <- introns.d[subjectHits(hits)]
bl.int <- blip[queryHits(hits)]
int.res <- cbind(as.data.frame(bl.int),as.data.frame(int.bl,row.names=c(1:length(int.bl))))
int.res$elem_ID <- paste(int.res[,1],int.res[,2],int.res[,3],sep="_") ##intronic blips
blip.b <- blip[-queryHits(hits)]

##Find nearest neighbor to the remaining blips
message("Identifying intergenic (non-intronic) novel elements")
nn <- distanceToNearest(blip.b,exon,ignore.strand=TRUE)
nn <- nn[nn@elementMetadata@listData$distance<10000]
nn.subj <- exon[subjectHits(nn)]
nn.qr <- blip.b[queryHits(nn)]
blip.c <- blip.b[-queryHits(nn)] ##blips far away from any RefSeq gene
ex.res <- cbind(as.data.frame(nn.qr),as.data.frame(nn.subj,row.names=1:length(nn.subj)))
ex.res$elem_ID <- paste(ex.res[,1],ex.res[,2],ex.res[,3],sep="_") ##blips near exons
message("Successfully identification of intronic, near-exon (<10 kb from exon), and orphan (>10kb from exon) elements")
write.table(ex.res,paste0(od,"/RE_NOVEL_nearexon.txt"))
write.table(int.res,paste0(od,"/RE_NOVEL_intron.txt"))
write.table(as.data.frame(blip.c),paste0(od,"/RE_NOVEL_orphan.txt"))
}
annotateRE(gtf.file)

/wsu/home/eh/eh97/eh9735/mott-ngs/Serono_MSE/SRE_discovery_v28/

extendExon <- function(gtf,od=out.directory,sample.names=sample.names.file,sampleloc=samples.location,reqdist=20){
	samples <- read.table(sample.names,as.is=T,sep="\t")[,1]
	ex.res <- read.table(paste0(od,"/RE_NOVEL_nearexon.txt"))
	nearex <- as(ex.res[,c(1:7,12:15)],"GRanges")
	int.res <- as(read.table(paste0(od,"/RE_NOVEL_intron.txt")),"GRanges")
	orphan.res <- as(read.table(paste0(od,"/RE_NOVEL_orphan.txt")),"GRanges")
	exon <- read.table(paste0(od,"/gtf_exons.bed"),as.is=T,sep="\t")
	colnames(exon) <- c("seqnames","start","end","grl_name","exon_name","strand")
	exon <- as(exon,"GRanges")
	nn <- distanceToNearest(nearex,exon,ignore.strand=TRUE)
	nn <- nn[nn@elementMetadata@listData$distance<reqdist]
	nn.subj <- exon[subjectHits(nn)]
	nn.qr <- nearex[queryHits(nn)]
	message(paste(length(nn.qr),"(non-intronic) novel elements are within",reqdist,"from an exon"))

	environment(extendCVG) <- environment()
	fin <- lapply(1:length(samples),extendCVG)
	pinex <- function(x){return(fin[[x]][[1]])}
	exlist <- lapply(1:length(samples),pinex)
	environment(getex) <- environment()
	excomb <- as(do.call(rbind,lapply(1:nrow(exlist[[1]]),getex)),"GRanges")

	##reduce excomb, and replace the annotations 
	c.ex <- reduce(excomb) ##130, before was 141
	hits= findOverlaps(c.ex,excomb,type="any",ignore.strand=TRUE)
	c.ex.p1 <- c.ex[queryHits(hits)]
	c.ex.p2 <- excomb[subjectHits(hits)]
	c.ex.p1$grl_name <- c.ex.p2$grl_name
	c.ex.p1$exon_name <- c.ex.p2$exon_name
	excomb <- c.ex.p1[!(duplicated(c.ex.p1))]

	##compile modified exons and original exons
	hits= findOverlaps(exon,excomb,type="within",ignore.strand=TRUE)
	modex <- c(exon[-unique(queryHits(hits))],excomb)  ##final exon set

	##combine blips from the sample
	pinbl <- function(x){return(fin[[x]][[2]])}
	bllist <- lapply(1:length(samples),pinbl)
	bllist <- reduce(do.call("c",bllist))
	blldiff <- setdiff(bllist,excomb,ignore.strand=TRUE)

	##replace the original exons with their extended exon
	hits= findOverlaps(nn.qr, nearex,type="within",ignore.strand=TRUE, minoverlap = 10L)
	nearfin <- nearex[-subjectHits(hits)] ##1996 ##compile the blldiff and the remaining nearexon blips
	bllfin <- c(nearfin,blldiff,ignore.mcols=TRUE) ##lost 19 blips
	##find the nearest exon to blips
	# nearest(bllfin,modex,ignore.strand=TRUE)
	modex[nearest(bllfin,modex,ignore.strand=TRUE)] ##annotation info
	exfin <- cbind(as.data.frame(bllfin,row.names=1:length(bllfin)),as.data.frame(modex[nearest(bllfin,modex,ignore.strand=TRUE)],row.names=1:length(bllfin))[,6:7])
	nexfin <- exfin[,c(1:3,6:7,5)] ##final near-exon set of blips
	nexfin[,5] <- "NOVEL_10KB_EXON"
	nexfin$ID <- paste(nexfin[,1],nexfin[,2],nexfin[,3],sep="_")
	colnames(nexfin) <- c("seqnames","start","end","gene_id","elem_anno","strand","elem_ID")

	##generate the bed file for the final exon set
	exfin <- as.data.frame(modex,row.names=1:length(modex))[,c(1:3,6:7,5)]
	exfin$ID <- paste(exfin[,1],exfin[,2],exfin[,3],sep="_")
	colnames(exfin) <- c("seqnames","start","end","gene_id","elem_anno","strand","elem_ID")

	##generate bed file for intron
	intfin <- as.data.frame(int.res)[,c(1:3,13,7,5)]
	intfin$ID <- paste(intfin[,1],intfin[,2],intfin[,3],sep="_")
	intfin[,5] <- "NOVEL_INTRONIC"
	colnames(intfin) <- c("seqnames","start","end","gene_id","elem_anno","strand","elem_ID")

	##generate bed file for orphans
	orfin <- data.frame(orphan.res)[,c(1:3,6:7,5)]
	orfin$ID <- paste(orfin[,1],orfin[,2],orfin[,3],sep="_")
	orfin[,5] <- "NOVEL_ORPHAN"
	orfin[,4] <- "NA"
	colnames(orfin) <- c("seqnames","start","end","gene_id","elem_anno","strand","elem_ID")

	p1 <- rbind(rbind(rbind(nexfin,intfin),orfin),exfin)
	write.table(p1,paste0(od,"/RE_compiled.txt"),row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
}
extendExon(gtf.file)

##As this section uses biomaRt, you must provide the biomart and dataset needed for the "useMart" function
annotateFinal <- function(extended,od=out.directory,mymart="ensembl",mydataset="hsapiens_gene_ensembl"){
	ensembl=useMart(biomart=mymart, dataset=mydataset)

	if(extended==TRUE){
		bd <- read.table(paste0(od,"/RE_compiled.txt"),as.is=TRUE,sep="\t")
		colnames(bd) <- c("seqnames","start","end","gene_id","class","strand","elem_ID")
		bd.tmp <- unique(unlist(lapply(strsplit(unique(unlist(strsplit(bd$gene_id,","))),"\\."),function(x) return(x[1])))) ##Get the total unique ensembl names
		# strsplit(unique(unlist(strsplit(bd$gene_id[1:5],","))),"\\.")
		bm <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","external_gene_name"),filters='ensembl_gene_id', values=bd.tmp,mart=ensembl)
		message("Removing suffix from ensembl gene id")
		bd$simpgene <- sapply(bd$gene_id,function(x) return(paste(unique(unlist(lapply(strsplit(unlist(strsplit(x,",")),"\\."),function(x) return(x[1])))),collapse=",")))
		message("Retreiving common gene name")
		message(paste("This step will take approximately",round(nrow(bd)/1000*5.717/60/60,digits=2),"hours.")) ##number of hours required to finish this step.
		bd$gene_name <- sapply(bd$simpgene,function(x) return(paste(bm[match(unlist(strsplit(x,",")),bm$ensembl_gene_id),"external_gene_name"],collapse=",")))
		fin <- bd[,c("seqnames","start","end","gene_name","class","strand","elem_ID","gene_id")]##reformat bd for bed file
		message("Saving the final elements as a bed file")
		write.table(fin,paste0(od,"/RE_complete.bed"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	}
	if(extended==FALSE){
		message("Formatting and combining your REs")
		ex.res <- read.table(paste0(od,"/RE_NOVEL_nearexon.txt"))
		nearex <- as(ex.res[,c(1:7,12:15)],"GRanges")
		int.res <- as(read.table(paste0(od,"/RE_NOVEL_intron.txt")),"GRanges")
		orphan.res <- as(read.table(paste0(od,"/RE_NOVEL_orphan.txt")),"GRanges")
		exon <- read.table(paste0(od,"/gtf_exons.bed"),as.is=T,sep="\t")
		colnames(exon) <- c("seqnames","start","end","gene_id","elem_anno","strand")
		exon$elem_ID <- paste(exon[,1],exon[,2],exon[,3],sep="_")

		##generate bed format for intron
		intfin <- as.data.frame(int.res)[,c(1:3,13,7,5)]
		intfin$ID <- paste(intfin[,1],intfin[,2],intfin[,3],sep="_")
		intfin[,5] <- "NOVEL_INTRONIC"
		colnames(intfin) <- c("seqnames","start","end","gene_id","elem_anno","strand","elem_ID")

		##generate bed format for orphans
		orfin <- data.frame(orphan.res)[,c(1:3,6:7,5)]
		orfin$ID <- paste(orfin[,1],orfin[,2],orfin[,3],sep="_")
		orfin[,5] <- "NOVEL_ORPHAN"
		orfin[,4] <- "NA"
		colnames(orfin) <- c("seqnames","start","end","gene_id","elem_anno","strand","elem_ID")

		##Generate bed format for near-exon
		nexfin <- data.frame(nearex)[,c(1:3,9,10,5)]
		nexfin$ID <- paste(nexfin[,1],nexfin[,2],nexfin[,3],sep="_")
		nexfin[,5] <- "NOVEL_10KB_EXON"
		colnames(nexfin) <- c("seqnames","start","end","gene_id","elem_anno","strand","elem_ID")

		bd <- rbind(exon,intfin,orfin,nexfin)
		bd.tmp <- unique(unlist(lapply(strsplit(unique(unlist(strsplit(bd$gene_id,","))),"\\."),function(x) return(x[1])))) ##Get the total unique ensembl names
		# strsplit(unique(unlist(strsplit(bd$gene_id[1:5],","))),"\\.")
		bm <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","external_gene_name"),filters='ensembl_gene_id', values=bd.tmp,mart=ensembl)
		message("Removing suffix from ensembl gene id")
		bd$simpgene <- sapply(bd$gene_id,function(x) return(paste(unique(unlist(lapply(strsplit(unlist(strsplit(x,",")),"\\."),function(x) return(x[1])))),collapse=",")))
		message("Retreiving common gene name")
		message(paste("This step will take approximately",round(nrow(bd)/1000*5.717/60/60,digits=2),"hours.")) ##number of hours required to finish this step.
		bd$gene_name <- sapply(bd$simpgene,function(x) return(paste(bm[match(unlist(strsplit(x,",")),bm$ensembl_gene_id),"external_gene_name"],collapse=",")))
		fin <- bd[,c("seqnames","start","end","gene_name","elem_anno","strand","elem_ID","gene_id")]##reformat bd for bed file
		message("Saving the final elements as a bed file")
		write.table(fin,paste0(od,"/RE_complete.bed"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	}
}
annotateFinal(TRUE)





