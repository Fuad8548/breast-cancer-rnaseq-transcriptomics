## Project 1 Results and Biological Interpretation

### Differential expression

A paired DESeq2 model (`~ patient_id + timepoint`) was applied to 48 RNA-seq samples from 24 matched patients with T1 and T3 measurements. After filtering, 13,903 genes were analyzed. The genome-wide T3-vs-T1 analysis identified 718 genes with nominal p < 0.05, 115 with p < 0.01, and 9 with p < 0.001, but no genes reached FDR < 0.05. The largest absolute estimated log2 fold change was approximately 2.09.

Therefore, individual genes were not classified as statistically significant DEGs on the basis of nominal p-values alone.

### Pathway-level analysis

Because coordinated changes can exist without individual genes surviving genome-wide multiple-testing correction, the complete DESeq2-ranked gene list was analyzed by preranked GSEA.

GO Biological Process GSEA revealed strong T3-associated enrichment for mitotic and genome-maintenance processes, including kinetochore/spindle function, chromosome segregation, cell-cycle progression, DNA replication, and homologous-recombination-related biology. T3-associated immune/lymphocyte terms were also enriched.

The T1 side showed enrichment involving cell/tissue organization and metabolic biology, including cell junction assembly, epithelial morphogenesis, extracellular/tissue structure, vascular development, carbohydrate metabolism, lipid metabolic regulation, and related processes.

Because GO terms are highly redundant, semantic simplification was applied before interpretation.

KEGG GSEA provided concordant pathway-level evidence. T3-associated pathways included Cell cycle, DNA replication, Homologous recombination, Ribosome, and Spliceosome. T1-associated pathways included glycan degradation and several signaling/tissue-interaction pathways such as integrin, focal-adhesion, gap-junction, cadherin, AMPK, and insulin signaling.

Disease-labelled KEGG pathways were interpreted as shared molecular signatures rather than evidence that the corresponding diseases were present in the samples.

### GO/KEGG convergence

Core/leading-edge genes from the strongest GO and KEGG pathways were intersected to identify genes repeatedly contributing to both annotation frameworks. This yielded 88 shared T3-associated core genes and 123 shared T1-associated core genes.

T3-associated enrichment was dominated by proliferative and genome-maintenance programs, including mitotic chromosome segregation, kinetochore and spindle-checkpoint activity, DNA replication, and homologous recombination. Core genes repeatedly contributing to these signals included KNL1, NDC80, SGO1, CDK1, CCNB1, BUB1/BUB1B, RAD51, BRCA2, and related genes.

T1-associated enrichment involved metabolic and tissue-interaction programs, including carbohydrate and small-molecule metabolism, integrin/focal-adhesion/cadherin signaling, gap-junction biology, and extracellular/tissue organization. Representative genes included GYS1, SCD, INSR, TNC, and related metabolic or extracellular/signaling genes.

### Patient-level candidate validation

Eight representative genes were examined across the 24 matched patients using paired VST expression. KNL1, GYS1, and INSR showed the clearest candidate-level paired evidence in the selected set, with 17 of 24 patients showing the expected directional shift and paired Wilcoxon p-values of 0.0135, 0.00894, and 0.0135, respectively. Other candidates showed the expected median direction but weaker within-patient consistency.

These candidate-level findings remain exploratory. They do not override the genome-wide DESeq2 conclusion that no individual gene passed FDR < 0.05, and they should not be described as validated biomarkers.

### Overall interpretation

Taken together, the analysis supports a model in which the T3 transcriptional state is characterized by coordinated enrichment of proliferative and genome-maintenance programs, especially mitotic chromosome segregation, kinetochore/spindle checkpoint activity, DNA replication, and DNA repair. The T1 state is comparatively characterized by metabolic and tissue-interaction-associated programs involving cell-matrix/cell-cell signaling, extracellular organization, and metabolism.

The analysis demonstrates an important distinction between gene-level and pathway-level inference: the dataset does not provide genome-wide FDR-significant individual DEGs under the fitted model, yet coordinated pathway-level shifts are detectable when the full ranked transcriptome is analyzed.
