# Exercice Pratique 1 : QC, Mapping et Variant Calling avec Conda

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

# Exercice Pratique 2 : QC, Mapping, Variant Calling et Automatisation avec Galaxy

## Objectif
Cet exercice a pour but de vous familiariser avec Galaxy pour réaliser un pipeline complet d’analyse NGS : **Contrôle Qualité (QC)**, **Trimming**, **Mapping (Alignement)**, et **Variant Calling**. Vous apprendrez également à automatiser ces étapes via un workflow pour optimiser les analyses répétées.

---

## Étape 1 : Configuration de l’Environnement

1. **Créer un Compte Galaxy** :
   - Rendez-vous sur [usegalaxy.org](https://usegalaxy.org) ou [usegalaxy.eu](https://usegalaxy.eu).
   - Cliquez sur **Register Here**, remplissez le formulaire et activez votre compte via l’email de confirmation.

2. **Créer un Historique** :
   - Connectez-vous et accédez au panneau d’historique à droite de l’interface.
   - Cliquez sur **+** pour créer un nouvel historique.
   - Renommez-le (ex. : `Exercice_Galaxy_QC_Mapping_VC`).

---

## Étape 2 : Importer les Données

1. **Importer les Données** :
   - Cliquez sur **Upload Data** dans le panneau de gauche.
   - Téléchargez les fichiers suivants :
     - Référence : `e_coli.fasta`.
     - Reads : `SRR2589044_1.fastq.gz` et `SRR2589044_2.fastq.gz`.

2. **Vérifier les Données** :
   - Assurez-vous que les fichiers sont bien visibles dans l’historique.
   - Vérifiez les formats (FASTA pour la référence, FASTQ pour les lectures).

---

## Étape 3 : Contrôle Qualité (QC)

1. **Exécuter FastQC** :
   - Recherchez **FastQC** dans la boîte à outils de gauche.
   - Sélectionnez les fichiers FASTQ des lectures.
   - Exécutez pour générer des rapports de qualité.

2. **Consolider les Rapports avec MultiQC** :
   - Recherchez **MultiQC**.
   - Sélectionnez les résultats FastQC.
   - Exécutez pour obtenir un rapport global.

   > **Résultat attendu :** Un rapport HTML consolidé des lectures brutes.

---

## Étape 4 : Nettoyage des Lectures

1. **Nettoyer les Lectures avec Fastp** :
   - Recherchez **Fastp** dans la boîte à outils.
   - Configurez les paramètres suivants :
     - Input : Fichiers FASTQ des lectures.
     - Qualité minimale : 20.
     - Longueur minimale des lectures : 25.
   - Exécutez l’outil.

2. **Valider le Nettoyage** :
   - Re-exécutez **FastQC** sur les fichiers nettoyés.
   - Compilez les rapports avec **MultiQC** pour confirmer les améliorations.

   > **Résultat attendu :** Un rapport MultiQC montrant une amélioration de la qualité des lectures après nettoyage.

---

## Étape 5 : Alignement (Mapping)

1. **Indexer la Référence** :
   - Recherchez **BWA-MEM Index** dans la boîte à outils.
   - Sélectionnez la référence `e_coli.fasta`.
   - Exécutez pour créer un index utilisable pour l’alignement.

2. **Aligner les Lectures avec BWA-MEM** :
   - Recherchez **BWA-MEM**.
   - Sélectionnez :
     - Référence indexée.
     - Lectures nettoyées (FASTQ appariés).
   - Exécutez pour générer un fichier SAM.

3. **Trier et Indexer les Alignements** :
   - Utilisez **SAMtools sort** pour convertir et trier le fichier SAM en BAM.
   - Exécutez **SAMtools index** pour indexer le fichier BAM.

   > **Résultat attendu :** Un fichier BAM trié et indexé, prêt pour l’appel de variantes.

---

## Étape 6 : Variant Calling

1. **Appeler les Variants avec FreeBayes** :
   - Recherchez **FreeBayes**.
   - Configurez les paramètres suivants :
     - Référence : `e_coli.fasta`.
     - Alignements : Fichier BAM trié et indexé.
   - Exécutez pour générer un fichier VCF brut.

2. **Filtrer les Variants avec BCFtools** :
   - Recherchez **BCFtools Filter**.
   - Appliquez les critères suivants :
     - Qualité (`QUAL`) supérieure à 20.
   - Exécutez pour produire un fichier VCF filtré.

   > **Résultat attendu :** Un fichier VCF contenant les variants de haute qualité.

---

## Étape 7 : Automatisation avec Workflow

1. **Créer un Workflow** :
   - Allez dans le menu **Workflow** et cliquez sur **Extract Workflow**.
   - Renommez et organisez les étapes pour inclure toutes les analyses effectuées (QC, trimming, mapping, variant calling).

2. **Exécuter le Workflow** :
   - Chargez un nouveau jeu de données dans un nouvel historique.
   - Exécutez le workflow pour automatiser toutes les étapes.

3. **Exporter le Workflow** :
   - Exportez le fichier `.ga` depuis Galaxy et incluez-le dans vos livrables.

---

## Questions

1. Quels problèmes de qualité avez-vous identifiés dans les lectures brutes ?
2. Comment les rapports FastQC/MultiQC ont-ils changé après le nettoyage ?
3. Quels types de variants (SNPs, indels) avez-vous détectés ?
4. Quels avantages offre l’utilisation d’un workflow automatisé ?

---

## Livrables

1. **Rapports MultiQC** :
   - Avant et après trimming.
2. **Fichier BAM** :
   - Trié (`sorted.bam`) et indexé (`sorted.bam.bai`).
3. **Fichier VCF** :
   - Brut (`raw.vcf`) et filtré (`filtered.vcf`).
4. **Workflow Galaxy** :
   - Fichier `.ga` exporté.
5. **Réponses aux Questions**.

