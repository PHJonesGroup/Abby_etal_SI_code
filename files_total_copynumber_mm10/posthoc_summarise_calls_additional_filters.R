# Apply some extra filters to remove artifacts

dat = readr::read_tsv("notch1_targeted_totalcn_full_segmentation.txt")
dat = dat[dat$chromosome!="X" & dat$width > 30000000 & dat$classification!="normal" & (dat$value > 0.1 | dat$value < -0.1),]
write.table(dat, file="notch1_targeted_totalcn_full_segmentation_filtered.txt", quote=F, sep="\t", row.names=F)
