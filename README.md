# HGT_experiment



# Coverage analysis of illumina samples

Scripts necessary to the analysis of HGT event between a pair of bacteria. 



1. Trimming
        # leading base quality 3
        # trailing base quality 3
        # phred quality of sequence, by sliding window of 4 bases and minimum phred score of 20
        # minimum length 50$
        # processing time <1m per paired reads.
       a=0;for i in $(ls *_R1_[0-9]*.fastq.gz | sort -t'_' -k3); do a=$((a + 1)); b=$(echo $i | cut -d'.' -f1 | cut -d'_' -f5); c=$(ls *_R2*$b*.fastq.gz); echo $i" "$c ;  java -jar /home/imateus/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $i $c $(echo $i | cut -d'_' -f1,2,3)_Trimm_Pair.fastq.gz $(echo $i | cut -d'_' -f1,2,3)_Trimm_Unpair.fastq.gz $(echo $c | cut -d'_' -f1,2,3)_Trimm_Pair.fastq.gz $(echo $c | cut -d'_' -f1,2,3)_Trimm_Unpair.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50;  done


        # ouput surviving paired reads: 
        ref VCA0107 Input Read Pairs: 2581550 Both Surviving: 2248355 (87.09%) Forward Only Surviving: 144415 (5.59%) Reverse Only Surviving: 61854 (2.40%) Dropped: 126926 (4.92%)
        ref VCA0747 Input Read Pairs: 2929340 Both Surviving: 2518201 (85.96%) Forward Only Surviving: 189556 (6.47%) Reverse Only Surviving: 73133 (2.50%) Dropped: 148450 (5.07%)

        # uncompress 
        # processing time < 5 sec
        for i in $(ls *Pair.fastq.gz); do echo $i ; gzip -d $i ;done
        
2. Mapping
- make index

           /home/imateus/software/novocraft/novoindex  -s 1 -n All_chromosomes-Sa5Y-VCA0107-frt-kan-frt All_chromosomes-Sa5Y-VCA0107-frt-kan-frt.nix All_chromosomes-Sa5Y-VCA0107-frt-kan-frt.fasta 

          /home/imateus/software/novocraft/novoindex  -s 1 -n All_chromosomes-Sa5Y-VCA0747-frt-kan-frt All_chromosomes-Sa5Y-VCA0747-frt-kan-frt.nix All_chromosomes-Sa5Y-VCA0747-frt-kan-frt.fasta



- mapping
                a=0;for i in $(ls *R1*Pair.fastq | sort -t'_' -k3); do a=$((a + 1)); b=$(echo $i | cut -d'.' -f1 | cut -d'_' -f5); c=$(ls *_R2*$b*.fastq); echo $i" "$c ; /home/imateus/software/novocraft/novoalign -d /home/imateus/Documents/V_cholerae/test1/ref_VCA0107/reference_template/*.nix -f $i $c -o SAM 2> 'stats'$i'VCA0107.txt' > $i'_VCA0107.sam'; done

                a=0;for i in $(ls *R1*Pair.fastq | sort -t'_' -k3); do a=$((a + 1)); b=$(echo $i | cut -d'.' -f1 | cut -d'_' -f5); c=$(ls *_R2*$b*.fastq); echo $i" "$c ; /home/imateus/software/novocraft/novoalign -d /home/imateus/Documents/V_cholerae/test1/ref_VCA0747/reference_template/*.nix -f $i $c -o SAM 2> 'stats'$i'VCA0747.txt' > $i'_VCA0747.sam'; done


- Filter to allowd only 3 mismatches only
        a=0;for i in $(ls *R1*Pair.fastq | sort -t'_' -k3); do a=$((a + 1)); b=$(echo $i | cut -d'.' -f1 | cut -d'_' -f5); c=$(ls *_R2*$b*.fastq); echo $i" "$c ;
        /home/imateus/software/novocraft/novoalign -d /home/imateus/Documents/V_cholerae/test1/ref_VCA0107/reference_template/*.nix -f $i $c -t 90 -o SAM 2> 'stats'$i'VCA0107.txt' > $i'_VCA0107.sam'; done



