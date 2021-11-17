.PHONY: clean

clean: 
	@python ./scripts/clean.py
	
detect: clean
	@python ./scripts/detect.py data/testing/scrambled.npy data/testing/mod_genome.fa data/testing/scrambled.for.bam Sc_chr04
	
reassemble: clean
	@python ./scripts/reassemble.py data/testing/scrambled.npy data/testing/mod_genome.fa data/testing/scrambled.for.bam Sc_chr04
	
train: 
	@python ./scripts/train.py

scramble: 
	@hiscram scramble test_data/seq.fa out_scrambles -1 test_data/sample.reads_for.fastq.gz -2 test_data/sample.reads_rev.fastq.gz --profile Debug
