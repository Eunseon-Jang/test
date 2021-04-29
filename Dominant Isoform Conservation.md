# Dominant Isoform Conservation

## 00)Tool/Raw file
### Tool

1. SRA toolkit : SRA data->fastq
2. flair-1.5 : Quantification
3. samtools : 용량이 큰 파일을 binary형태인 bam포맷으로 바꿈
4. bedtools : BED, GFF3, VCF등 1차원 좌표계를 포함하는 유전체 자료를 빠르게 상호계산 하도록함.
5. gffread : gtf -> gff
6. requirements : intervaltree, Cython, kerneltree(Cython이 있어야 깔수있음), tqdm, pybedtools, pysam, numpy
<br/>
root 권한이 필요한것들.(sudo) --user 옵션으로 깔기
<br/>
pip install intervaltree --user
<br/>
pip install kerneltree
<br/>
pip install Cython --user
<br/>
pip install kerneltree --user
<br/>
pip install tqdm --user
<br/>
pip install pybedtools --user
<br/>
pip install pysam --user
<br/>
pip install numpy --user
<br/>

### Raw file
~~1. iGenomes
https://support.illumina.com/sequencing/sequencing_software/igenome.html~~
<br/>
2. NCBI(genome,GFF)
<br/>
Homo sapiens : https://www.ncbi.nlm.nih.gov/genome/?term=Homo+sapiens
<br/>
Mus musculus : https://www.ncbi.nlm.nih.gov/genome/?term=mus+musculus

### 00.1)UCSC는 NM.xxxx만 있는데 NCBI는 NM,XR annotation둘다 있어서 바꿨음. NCBI는 gff 파일이므로 gft변경함
gffread -T -o GCF_000001405.39_GRCh38.p13_genomic.gtf GCF_000001405.39_GRCh38.p13_genomic.gff
<br/>
gffread -T -o GCF_000001635.27_GRCm39_genomic.gtf GCF_000001635.27_GRCm39_genomic.gff
<br/>

## 01)Alignment
~~python /home/nelljk98/tools/flair-1.5/flair.py align --genome /home/nelljk98/Dominant_isoform_conservation/00_Data_liver/Genome_data/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --reads /home/nelljk98/Dominant_isoform_conservation/00_Data_liver/Fastq/Mus_musculus/ERR2844020.1.fastq --threads 8 --output mouse_liver_FLAIR_align~~
<br/>
Mouse
<br/>
python /home/nelljk98/tools/flair-1.5/flair.py align --output mouse_liver_FLAIR_align --threads 8 --genome /home/nelljk98/Dominant_isoform_conservation/00_Data_liver/Genome_data/NCBI/Mus_musculus/GCF_000001635.27_GRCm39_genomic.fna --reads /home/nelljk98/Dominant_isoform_conservation/00_Data_liver/Fastq/Mus_musculus/ERR2844020.1.fastq 
<br/>
samtools sort -@ 8 -o mouse_liver_FLAIR_align.sorted.bam mouse_liver_FLAIR_align.bam
<br/>
bamToBed -bed12 -i mouse_liver_FLAIR_align.sorted.bam > mouse_liver_FLAIR_align.bed12
<br/>
Human
<br/>
samtools sort -@ 8 -o human_liver_FLAIR_align.sorted.bam human_liver_FLAIR_align.bam
<br/>
bamToBed -bed12 -i human_liver_FLAIR_align.sorted.bam > human_liver_FLAIR_align.bed12

## 02)Correct
python /home/nelljk98/tools/flair-1.5/flair.py correct --query ../../01_FLAIR/NCBI/mouse_liver_FLAIR_align.bed12 --genome /home/nelljk98/Dominant_isoform_conservation/00_Data_liver/Genome_data/NCBI/Mus_musculus/GCF_000001635.27_GRCm39_genomic.fna --gtf /home/nelljk98/Dominant_isoform_conservation/00_Data_liver/Genome_data/NCBI/Mus_musculus/GCF_000001635.27_GRCm39_genomic.gtf --output mouse_liver_FLAIR_correct --threads 8

## 03)Collapse
python /home/nelljk98/tools/flair-1.5/flair.py collapse --query ../02_FLAIR_correct/NCBI/mouse_liver_FLAIR_correct_all_corrected.bed --genome /home/nelljk98/Dominant_isoform_conservation/00_Data_liver/Genome_data/NCBI/Mus_musculus/GCF_000001635.27_GRCm39_genomic.fna --gtf /home/nelljk98/Dominant_isoform_conservation/00_Data_liver/Genome_data/NCBI/Mus_musculus/GCF_000001635.27_GRCm39_genomic.gtf --reads /home/nelljk98/Dominant_isoform_conservation/00_Data_liver/Fastq/Mus_musculus/ERR2844020.1.fastq --support 1 --isoformtss --trust_ends --keep_intermediate --generate_map --output mouse_liver_FLAIR_collapse
## 04)Quantify
python /home/nelljk98/tools/flair-1.5/flair.py quantify --reads_manifest reads_manifest.mouse.tsv --isoforms ../03_FLAIR_collapse/mouse_liver_FLAIR_collapse.isoforms.fa --threads 8 --tpm --trust_ends --output mouse_liver_FLAIR_quantify
## 05)Edit file
CompleteAnnotation, Partial Annotation, Novel Annotation 파일로 각각 나눔 

**Transcript ID, Gene ID 모두 가진 Isoform(해당 isoform이 알려진) (Complete Annotation)**. 

cat human_brain_FLAIR_quantify.tpm.tsv |grep "gene-" |grep "rna-" > human_brain_FLAIR_quantify.tpm.CompleteAnnotation.txt  

**Transcript ID는 없지만 Gene ID를 가진 Isoform(어떤 유전자의 isoform인지는 알지만 알려지지 않은 isoform) (Partial Annotation)**. 

cat human_brain_FLAIR_quantify.tpm.tsv |grep "gene-" |grep -v "rna-" > human_brain_FLAIR_quantify.tpm.PartialAnnotation.txt 

**특정 유전자 위치의 Isoform이 아님 (Novel Annotation)**. 

cat human_brain_FLAIR_quantify.tpm.tsv |grep -v "gene-" |grep -v "rna-" > human_brain_FLAIR_quantify.tpm.NovelAnnotation.txt
