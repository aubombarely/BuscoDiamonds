# BuscoDiamonds

  This program has several functions:

   1. Evaluate the completeness of the gene space in a whole genome sequencing dataset that has not been assembled.
   2. Reconstruct the gene/CDS/protein associated with the hits of the Gene Space reference.
   3. Evaluate the duplication level of the analyzed Gene Space.
   4. Estimate the genome coverage and heterozygosity based in the reconstruction of the evaluation Gene Space set.

  We define the gene space as the gene population in the genome of a specific
  individual or accession. So, for example, the gene space for the Arabidopsis
  thaliana Col. genome (annotation TAIR10) is composed by 27,416 genes.

  A gene space representative is serie of proteins that could be found in most
  of the species of a specific clade. So for example, if all the Brassicaceae
  species shared the same 8,000 genes, a gene space representaive for the
  Brassicaeae clade will have an ancestor single copy ortholog of those genes.
  You can find more information about the ortholog datasets at Zdobnov et al.
  2017 (OrthoDB v9.1) and Waterhouse et al. 2017 (BUSCO3).  
  
  This scripts uses DIAMOND (Buchfink et al. 2015) to align unassembled short
  reads against a BUSCO/OrthoDB dataset to evaluate the completness of the
  genome before the assembly, but at the same time re-construct the genes, 
  CDS and proteins to retrieve more information (e.g. heterozygosity, coverage
  ...).

  The script uses two mandatory parameters:

   1. Reference protein sequence file (-r)
   2. A list of fastq short read files to map to the reference (-i). Multiple files can be specified using the comma. The script DOES NOT USE the pair end information

  Prerrequisites
  ==============
   * Bioperl
   * Diamond accessible from the path or as the env. variable DIAMOND_PATH
   * Seqtk accessible from the path or as the env. variable SEQTK_PATH
   * Minia accessible from the path or as the env. variable MINIA_PATH
   * Bowtie2 accessible from the path or as the env. variable BOWTIE2_PATH
   * Mapsembler2 accessible from the path or as the env. variable MAPSEMBLER_PATH
   * Augustus accessible from the path or as the env. variable AUGUSTUS_PATH
     Additionally the the Augustus config directory should be accessible as the
     environmental variable AUGUSTUS_CONFIG_PATH
   * Samtools accessible from the path or as the env. variable SAMTOOLS_PATH
   * Bedtools accessible from the path or as the env. variable BEDTOOLS_PATH 
   * Freebayes accessible from the path or as the env. variable FREEBAYES_PATH

  Steps
  =====
  ```
     +---- 
     | 1.1- Create the DIAMOND database.
     | 1.2- Align the reads against the DIAMOND database.
   1x| 1.3- Filter the DIAMOND output if a filter was specified
     | 1.4- Estimate the completeness of the REF_GENE_SPACE based in the mapping
     | 1.5- Select of the Fastq reads that were mapped to the reference
     +-----
     +-----
     | 2.1- Run Minia assembler on the Fastq reads
   Rx| 2.2- Extend the Minia contigs with Mapsembler2
     | 2.3- Mapping the whole read dataset back to the extended contigs
     | 2.4- Extract the reads from the mapping
     +-----
   1x| 2.5- Reassemble the fastq reads with Minia
     +-----  
     +-----
     | 3.1- Split contigs in separate files
   1x| 3.2- Run Augustus in a multithread mode
     | 3.3- Parse the Augustus output to get CDS and protein sequences
     +-----
     +-----
     | 4.1- Evaluate the completness in the assembly (final contigs) with BUSCO
     | 4.2- Evaluate the completness in the predicted proteins with BUSCO
   1x| 4.3- Reads remapping with the assembly
     | 4.4- Remapping file sorting and indexing
     | 4.5- Coverage analysis using Bedtools
     | 4.6- Variant calling usinf Freebayes to evaluate heterozygosity
     +-----
     Note1: The number of assembly rounds (R) will be defined by the parameter -a <assembly rounds>
  ```
   
  Arguments for External Tools
  ============================
  BuscoDiamonds can pass some arguments to some tools such as:
   * Minia: -m <minia_arguments> (e.g. -m '-abundance-min=5')
   * Mapsembler2: -e <extremities:args;extend:args> (e.g. -e 'extend:-c=5')
   * Bowtie2: -w <bowtie2_arguments> (e.g. -w '--very-fast')

   - Note1: Threading options will be overwritten by -t <threads>
   - Note2: Imput/Output/Formating options are banned.
   - Note3: To use Kmers > 31 for Minia and Mapsembler2, you may need to recompile the tool using make k=64

  Completeness Evaluation
  =======================
  BuscoDiamonds evaluates the completeness of the reference Gene Space using two
  different approaches:

     A- Analyzing the reads that align with the protein model (e.g. BUSCO protein)
        and calculating if they cover the length (+/- SD) of the protein model
        at the step 1.5. 
        
        Because the assembly has not been performed yet, it will only estimate
        if the protein is represented in a "complete" or in a "fragmented" form.
        Option -l <std_length_file> controls the SD values for the protein models.
        If no option -l is supplied, it will use 10% of the protein length. No 
        expected score is used, but the Diamond hits can be filtered prior this 
        analysis with the -f option (e.g. -f 'ZI=60' to keep 60% of similarity 
        or more).

     B- BUSCO run over the assembled contigs (step 2.5) using the BUSCO mode 
        'genome' and over the predicted proteins (step 3.4) using the 'protein'
        mode. This option is only available with -b <BUSCO_lineage> option is 
        supplied.

  Sequence Coverage Evaluation
  ============================
  Sequence coverage evaluation is performed mapped the reads back to the contigs
  and running Bedtools genomecov on the bam files. The script will calculate
  the minimum, maximum and mean for each contig and for the whole assembly 
  (marked as WholeAssembly) and report them in the FinalAssembly.coverage.txt 
  file

  Heterozygosity
  ==============
  BuscoDiamonds will also evaluate the heterozygosity of the sample calling
  variants on the bam file using Freebayes.

  Output Directories and Files
  ==========================
  -------------------------------------------------------------------------
  00_DiamondReference: Directory with the Diamond database (*.dmnd file)
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  01_DiamondMapping: Directory with the results of the Diamond mapping
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    DiamondMapping.CompSummary.txt   
     - Gene space evaluation with seven columns: GroupID, ExpectedLength, 
       ExpectedLengthSD (preset by -l option or the 10% of the ExpectedLength),
       NumberOfReadHits, AaCoveredByOneOrMoreReads, AaCoveredByFiveOrMoreReads 
       and Tag classifying the coverage of the RefSeq.    
    DiamondMapping.coverage.txt      
     - Mapping coverage file with one line per reference sequence. Each line has
       three columns: RefSeqID, RefSeqLength, MappingDepthByPosition 
       (separated by commas) 
    DiamondMapping.PepRefSummary.txt 
     - Gene space evaluation summary based on the Diamond mapping
    DiamondMapping.sam               
     - Mapping file in Sam format (Diamond output)
  -------------------------------------------------------------------------
  02_ExtractedReads: Directory for the extracted reads from the Diamond mapping
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    DndHitsReadIDs.txt
     - List of the read IDs from the Diamond output
    DndHitsReadSet.fq
     - Reads mapped to the Reference File
  -------------------------------------------------------------------------
  03_Assembly: Directory with one subdirectory per round (Assembly_RoundX)
               where X is the assembly round. It also will contain the final
               assembly file (MergedMiniaAssembly.contigs.fa), the formated
               assembly file with contig_ids as CONTIGXXXXXXX where X is the
               contig number (FinalAssembly.contigs.fa) and the indexes of the
               final assembly for the Bowtie2 mapping. The Assembly_RoundX
               directory contains:
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MergedMiniaAssembly.contigs.fa
     - Contigs for Minia assembly
    MergedMiniaAssembly.contigs.h5
     - Kmer graphs for Minia assembly
    MergedMapsember2Extremities.fasta
     - Extremities for the Mapssembler2 run
    MergedMapsember2Extend_k_K_c_C_t_2.fasta
     - Extended contigs for Mapssembler2 run where K=Kmer, C=Coverage and 
       This file will be renamed to the file MergedExtendedContigs.fa
    MergedExtendedContigs.*.bt2
     - Indexes for Bowtie2
    MergedMappedReads.sam
     - Mapped reads to the MergedMapsember2Extend_k_K_c_C_t_2.fasta file
    MergedExtractedReads.fastq
     - Extracted reads from the MergedMappedRead.sam file  
  -------------------------------------------------------------------------
  04_StructuralAnnotation: Directory with the structural annotation of the 
                           contigs. It will have three subdirectories
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    01_OneContigPerFile
     - A directory with single sequence files to multithread Augustus
    02_AugustusPredictions
     - A directory with the one file per Augustus prediction.
    03_GeneModels
     - A directory with three files: 1- MergedAugustus.gff with the gene models
       2- MergedAugustus.CDS.fasta, fasta file with the predicted CDS and 3-
       MergedAugustus.PEP.fasta, a fasta file with the predicted proteins 
 -------------------------------------------------------------------------
  05_GeneSpaceEvaluation: Directory with three subdirectories with the BUSCO
                          run for the contigs (01_GeneSpace_By_Genome) and 
                          the proteins (02_GeneSpace_By_ProteinSet). It also
                          contains the 03_ReadRemappingEvaluation directory
                          with the following files:
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    FinalCtgsMappedReads.sam/bam/sorted.bam/sorted.bam.bai
      - SAM/BAM files after the mapping of the whole read dataset to the contigs
    FinalAssembly.coverage.txt
      - Bed file with the coverage per position (output of Bedtools genomecov)
    FinalAssembly.SelfmapVariants.vcf
      - VCF file with the variant calling (Freebayes on the remapped BAM)
    FinalAssembly.GFFvsVCF.bed
      - GFF/VCF intersect Bed file  (output of Bedtools insersect)
    SummaryGenomeCoverage.txt
      - Summary with the sequencing coverage analysis (with 4 columns: 
        ContigID, minimum coverage, Maximum coverage and mean coverage). 
        The WholeAssembly is printed in the second line.
    SummaryGenomeVariants.txt
      - Summary with the variants, with the following columns: SeqID/FeatID, 
        type, length, mean coverage/depth, heterozygous variants, homozygous 
        variants, SNPs, insertions, deletions, mnp and complex variants.
  -------------------------------------------------------------------------
  06_Finalresults: Directory with the final results
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    BuscoDMD_CDS.fasta
     - Predicted CDS with the SeqID as Prefix (-p option) + RefSeqIDMatch + N 
       (for multiple copies). In the description, it is presented the original 
       geneID assigned by Augustus.
    BuscoDMD_PEP.fasta
     - Predicted proteins (with the same SeqIDs than the CDS file)
    BuscoDMD_RunSummary.txt
     - Summary of the run with stats and files.
  -------------------------------------------------------------------------


  Examples of running commands
  ============================
  Basic command without run BUSCO
  ``` 
   BuscoDiamonds -r BuscoPlantProteins.fasta -i MySp_R1.fq,MySp_R2.fq -o MySp_BuscoDiamonds2018 -l BuscoPlantProtein_lengths_cutoff
  ```

  Basic command running BUSCO
  ```  
    BuscoDiamonds -r BuscoPlantProteins.fasta -i MySp_R1.fq,MySp_R2.fq -o MySp_BuscoDiamonds2018 -l BuscoPlantProtein_lengths_cutoff -b BuscoDirectoryWithTheLineage
  ```

  Filtering the hits with percentage of similarity below 60%
  ```
    BuscoDiamonds -r BuscoPlantProteins.fasta -i MySp_R1.fq,MySp_R2.fq -o MySp_BuscoDiamonds2018 -f ZI=60
  ```
