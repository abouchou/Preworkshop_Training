Je rencontre des difficultés pour générer un lien de téléchargement pour le fichier Markdown. Je vais réessayer une méthode différente pour vous fournir le fichier. Un instant, s'il vous plaît.

Je ne parviens pas à générer un fichier téléchargeable pour le Markdown. Je vous propose de copier directement le contenu suivant dans un fichier `.md` :

---

```markdown
# Exercice Pratique : QC, Mapping et Variant Calling avec Conda

## Objectif
L’objectif de cet exercice est de vous permettre d’appliquer les concepts appris en bioinformatique pour exécuter les étapes clés d’un pipeline d’analyse NGS dans un environnement Conda. Vous utiliserez des outils comme FastQC, Trimmomatic, BWA, SAMtools, et FreeBayes.

## 1. Installation de l’Environnement

1. Installez Miniconda si ce n'est pas encore fait, en suivant [ce lien](https://docs.anaconda.com/miniconda/).
2. Créez un environnement nommé `analyse` :
    ```bash
    conda create -n analyse -y
    conda activate analyse
    ```
3. Ajoutez les channels nécessaires et installez les outils :
    ```bash
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda install fastqc multiqc trimmomatic bwa samtools freebayes bcftools -y
    ```

---

## 2. Téléchargement des Données

1. Créez une structure de répertoires :
    ```bash
    mkdir -p ~/analyse/{reference,raws,results}
    cd ~/analyse
    ```
2. Téléchargez les fichiers suivants :
   - Référence (E. coli) :
        ```bash
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz -O reference/e_coli.fasta.gz
        gunzip reference/e_coli.fasta.gz
        ```
   - Reads (FASTQ) :
        ```bash
        cd raws
        wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
        wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
        wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/001/SRR2584861/SRR2584861_1.fastq.gz
        wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/001/SRR2584861/SRR2584861_2.fastq.gz
        ```

---

## 3. Contrôle Qualité (QC)

1. Exécutez FastQC sur les fichiers bruts :
    ```bash
    mkdir ../results/fastqc_raw
    fastqc *.fastq.gz -o ../results/fastqc_raw/
    ```
2. Compilez les rapports avec MultiQC :
    ```bash
    cd ../results/fastqc_raw
    multiqc .
    ```

---

## 4. Trimming des Lectures

1. Créez un fichier d’adaptateurs `NexteraPE-PE.fa` avec le contenu suivant :
    ```
    >PrefixNX/1
    AGATGTGTATAAGAGACAG
    >PrefixNX/2
    AGATGTGTATAAGAGACAG
    ```
2. Lancez Trimmomatic pour nettoyer les lectures :
    ```bash
    cd ../../raws
    mkdir ../results/trimming
    trimmomatic PE SRR2589044_1.fastq.gz SRR2589044_2.fastq.gz \
    ../results/trimming/SRR2589044_1_paired.fq.gz ../results/trimming/SRR2589044_1_unpaired.fq.gz \
    ../results/trimming/SRR2589044_2_paired.fq.gz ../results/trimming/SRR2589044_2_unpaired.fq.gz \
    ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:25
    ```

---

## 5. Alignement (Mapping)

1. Indexez la Référence :
    ```bash
    cd ../reference
    bwa index e_coli.fasta
    ```
2. Effectuez le Mapping avec BWA-MEM :
    ```bash
    mkdir ../results/mapping
    bwa mem -R "@RG\\tID:SRR2589044\\tSM:SRR2589044\\tPL:ILLUMINA" e_coli.fasta \
    ../results/trimming/SRR2589044_1_paired.fq.gz ../results/trimming/SRR2589044_2_paired.fq.gz > ../results/mapping/SRR2589044.sam
    ```

---

## 6. Conversion et Tri des Fichiers

1. Convertissez le fichier SAM en BAM et triez-le :
    ```bash
    cd ../results/mapping
    samtools view -S -b SRR2589044.sam > SRR2589044.bam
    samtools sort SRR2589044.bam -o SRR2589044_sorted.bam
    samtools index SRR2589044_sorted.bam
    ```

---

## 7. Variant Calling

1. Appelez les variants avec FreeBayes :
    ```bash
    mkdir ../variant_calling
    freebayes -f ../../reference/e_coli.fasta SRR2589044_sorted.bam > ../variant_calling/SRR2589044.vcf
    ```
2. Filtrez les variants avec BCFtools :
    ```bash
    bcftools filter -s LowQual -e 'QUAL<20 || DP<10' ../variant_calling/SRR2589044.vcf > ../variant_calling/SRR2589044_filtered.vcf
    ```

---

## Questions

1. Quels types de variants (SNPs, Indels) avez-vous observés ?
2. Comparez les fichiers VCF avant et après filtrage : quelles différences remarquez-vous ?
3. Pourquoi est-il important de vérifier la qualité des données avec FastQC/MultiQC avant le mapping ?

---

## Livrables

- Les fichiers suivants doivent être soumis :
  1. `multiqc_report.html` pour le QC.
  2. `SRR2589044_sorted.bam` pour le mapping.
  3. `SRR2589044_filtered.vcf` pour l’appel de variantes.
```

---

