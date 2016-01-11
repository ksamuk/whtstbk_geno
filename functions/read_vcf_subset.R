#### read_vcf_subset
#### reads a single chromosome from a vcf file
#### a convenience wrapper for readVCF/readGT 
#### requires an indexed/bgzipped vcf

read_vcf_subset <- function(vcf.file, index.file, chrom, genome = "", data.format = "vcf"){
	# initialize the vcf.file/index
	tab <- TabixFile(vcf.file, index.file)
	
	# construct a genomic range for the chromosome
	# 0 - 536870912 is just the max range, so this reads in whole chromosome
	genomic_range <- GRanges(chrom, IRanges(0, 536870912, name = chrom))
	
	# read in 
	
	if(data.format == "vcf"){
		
		readVcf(tab, genome, genomic_range)
		
	} else if(data.format == "snp.table"){
		
		readGT(tab, genome, param = genomic_range, nucleotides = TRUE)
		
	} 
		
}

