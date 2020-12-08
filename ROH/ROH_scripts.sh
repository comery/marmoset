#Step 1: Reads mapping and alignment filtering
#bwa-mem
bwa mem -t 6  $ref in_1.clean.fq.gz in_2.clean.fq.gz | samtools view -q 10 -Sb - >  out.bam 

#filtering
samtools sort -@ 6 -m 2G -O bam -o sorted.bam  out.bam 
samtools index sorted.bam 
samtools view  -h  -F 12 -F 256 sorted.bam |grep -v 'XA:Z' | grep -v 'SA:Z' |$samtools view -b - >sorted.uniq.bam
samtools addreplacerg  -r '@RG\tID:in\tPL:illumina\tSM:in' -o sorted.uniq.addRG.bam sorted.uniq.bam
samtools index sorted.uniq.addRG.bam

#Step 2: Genome vcf calling
#genome vcf calling
$gatk HaplotypeCaller \
        -R $ref \
        --emit-ref-confidence GVCF \
        -I sorted.uniq.addRG.bam \
        --native-pair-hmm-threads 10 \
        -O out.g.vcf 


#Step 3-1: SNP calling for each individual for Depth calculation
#GenotypeGVCF
$gatk GenotypeGVCFs \
        -R $ref \
        -V out.g.vcf \
        -O out.vcf 
# filter snp
$gatk SelectVariants \
                -select-type SNP \
                -V out.vcf \
                -O out.snp.vcf.gz
$gatk VariantFiltration \
                -V $a.snp.vcf.gz \
                --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
                --filter-name "Filter" \
                -O out.snp.filter.vcf.gz
##draw plot for each individual after obtain the depth
library(ggplot2)
a<-read.table("input.txt",sep="\t",header=F)
b<-log(a$Depth)
c<-cbind(a,b)
colnames(c)<-c("Sample","Chromosome","Location","Depth","LogNum")
pal<-c('#E9B7ED','#f38181','#fce38a','#00b8a9')
c$Chromosome<- factor(c$Chromosome,levels=c("Chr2","Chr3","Chr4","Chr13"))
ggplot(data=c,aes(x=Location,y=LogNum,fill=Sample,color=Sample,group=Sample))+geom_line()+facet_grid(Chromosome~.)+xlab("Location")+ylab("Log(Variants_Number)")+scale_size_manual(values=c(0.1, 0.1,0.1,0.1))+ scale_colour_manual(values = pal)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_set(theme_bw())
dev.off()


#Step 3-2: Roh analysis using all g.vcf
#combine gvcfs for 11 individuals
list_of_gvcfs="gvcf_11samples.lis"
##"gvcf_11samples.lis": the list file with 11 IDs.out.g.vcf
gvcf_argument=""
while read -r line; do
gvcf_argument=$gvcf_argument" -v $line"
done < "$list_of_gvcfs"

$sentieon driver -r $ref   --algo GVCFtyper  \
 --call_conf  30 --emit_conf 30  \
 --emit_mode variant  \
 --max_alt_alleles 2 \
 --annotation 'AC,AF,AN,BaseQRankSum,QD,FS,SOR,MQ,MQRankSum,ReadPosRankSum' \
 $gvcf_argument joint.vcf.gz

#SNP calling for 11 samples
time gatk-4.0.5.2/gatk SelectVariants \
    -select-type SNP \
    -V joint.vcf.gz \
    -O snp.vcf.gz
time gatk-4.0.5.2/gatk SelectVariants \
    -select-type INDEL \
    -V joint.vcf.gz \
    -O indel.vcf.gz
time gatk-4.0.5.2/gatk VariantFiltration \
    -V indel.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "HARD_fail" \
    -O indel.filter.vcf.gz

time gatk-4.0.5.2/gatk VariantFiltration \
    -V snp.vcf.gz \
    --filter-expression "AN < 22" --filter-name "AN_FAIL" \
    --filter-expression "QUAL < 30.0" --filter-name "QUAL_FAIL" \
    --filter-expression "DP < 50.0 || DP > 1500.0" --filter-name "DP_FAIL" \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 4.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"  --filter-name "HARD_FAIL" \
    --mask indel.filter.vcf.gz --mask-extension 3 --mask-name aroundIndel \
    --cluster-size 3 --cluster-window-size 10 \
    -O snp.filter.vcf.gz
time gatk-4.0.5.2/gatk SelectVariants \
    -V snp.filter.vcf.gz \
    --exclude-filtered \
    -O snp.pass.vcf.gz

vcftools --max-missing 1  --recode-INFO-all  --recode --out final --gzvcf snp.pass.vcf.gz   --min-alleles 2 --max-alleles 2  --remove-filtered-all   --stdout | bgzip -c > final.snp.vcf.gz
vcftools --FILTER-summary  --out check --gzvcf final.snp.vcf.gz

#combine all 12 samples
bcftools-1.8/bin/bcftools merge --output combine.vcf.gz --output-type z --missing-to-ref  --threads 2 final.snp.vcf.gz  offspring.vcf.gz
#keep only normal chromosomes
gzip -dc combine.vcf.gz  |grep -v ChrX | grep -v scaffold |bgzip -c >combine_noX.vcf.gz
#maf & hwe
vcftools --max-missing 1  --maf 0.2  --hwe 0.001  --recode-INFO-all  --recode  --gzvcf  combine_noX.vcf.gz  --stdout |bgzip -c > use.vcf.gz
#ROH
/hwfssz1/ST_DIVERSITY/PUB/USER/zhouchengran/software/bcftools-1.8/bcftools roh  -O r --AF-dflt 0.4 -o roh  use.vcf.gz
awk '$7>=10' roh>roh.fil.results

