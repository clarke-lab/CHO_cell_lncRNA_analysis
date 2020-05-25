gffread lncrna_annotation/lncrna.gtf -g reference_genome/ensembl_chok1_genome.fa  -w plot_data/chok1.lncrna.fasta?
gffread transcriptome_assembly/stringtie.protein.coding.gtf -g reference_genome/ensembl_chok1_genome.fa  -w plot_data/chok1.mrna.fasta?
gffread transcriptome_assembly/stringtie.all.transcripts.gtf -g reference_genome/ensembl_chok1_genome.fa  -w plot_data/chok1.transcripts.fasta?
bioawk -c 'fastx' '{print gc($seq)*100"\t""lncRNA"}' plot_data/chok1.lncrna.fasta > plot_data/gc.content.txt
bioawk -c 'fastx' '{print gc($seq)*100"\t""mRNA"}' plot_data/chok1.mrna.fasta > plot_data/gc.content.txt
 

bioawk -c 'fastx' '{print gc($seq)*100"\t""lncRNA"}' plot_data/chok1.lncrna.fasta? >


ioawk -c 'fastx' '{print gc($seq)*100"\t""lncRNA"}' plot_data/chok1.lncrna.fasta? >> mgc.content.txt
bioawk -c 'fastx' '{print gc($seq)*100"\t""mRNA"}' plot_data/chok1.mrna.fasta? >> mlot_data/igure2/gc.content.txt




