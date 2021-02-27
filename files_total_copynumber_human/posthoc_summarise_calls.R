# Here we combine all the output profiles from the totalcn analysis (from a listing to be put in segment_files.lst)

library(GenomicRanges)

infiles = read.table("segment_files.lst", header=F, stringsAsFactors=F)$V1
infiles = infiles[!grepl("output_old_pipeline", infiles)]
dat = lapply(infiles, function(x) { sam=unlist(strsplit(basename(x), "_"))[1]; temp=read.table(x, header=T, stringsAsFactors=F); temp$samplename=sam; temp; })
dat = do.call(rbind, dat)
dat = makeGRangesFromDataFrame(dat, keep.extra.columns=T)
notch1 = makeGRangesFromDataFrame(data.frame(chromosome=2, start=26457903, end=26516663))
overlap = findOverlaps(notch1, dat)
dat_notch1 = dat[subjectHits(overlap),]

cast_to_df = function(dat) {
dat = as.data.frame(dat)
colnames(dat)[1] = "chromosome"
dat = dat[,!grepl("strand", colnames(dat))]
return(dat)
}

dat = cast_to_df(dat)
dat_notch1 = cast_to_df(dat_notch1)
# make a complete table, even when no data is available for notch1
mapping = read.table("../sample_mapping.lst", header=F, stringsAsFactors=F)
missing = mapping$V1[!mapping$V1 %in% dat_notch1$samplename]
missing = missing[missing %in% dat$samplename]
for (samplename in missing) {

segment_before = tail(dat[dat$samplename==samplename & dat$chromosome=="2" & dat$start < 26457903,], 1)
segment_after = head(dat[dat$samplename==samplename & dat$chromosome=="2" & dat$end > 26516663,], 1)

temp = dat[1,,drop=F]
temp[,1:ncol(temp)] = NA
temp$samplename = samplename
temp$classification = "no data"
if (segment_before$classification=="loss") {
	temp$classification = paste0(temp$classification, " - loss before")
}
if (segment_before$classification=="gain") {
	temp$classification = paste0(temp$classification, " - gain before")
}
if (segment_after$classification=="loss") {
	temp$classification = paste0(temp$classification, " - loss after")
}
if (segment_after$classification=="gain") {
	temp$classification = paste0(temp$classification, " - gain after")
}

dat_notch1 = rbind(dat_notch1, temp)
}

write.table(dat, file="notch1_targeted_totalcn_full_segmentation.txt", quote=F, sep="\t", row.names=F)
write.table(dat_notch1, file="notch1_targeted_totalcn_onlynotch1_segmentation.txt", quote=F, sep="\t", row.names=F)
