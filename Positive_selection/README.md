### identification of orthologous genes

We used BLAST Reciprocal Best Hits (RBH) method (Supplementary Note) to identify high-confidence one-to-one orthologous genes among species, including other three NWM (white-faced capuchin (Cebus capucinus), Ma's night monkey (Aotus nancymaae), black-capped squirrel monkey (Saimiri boliviensis)) and three old world primates (human (Homo sapiens), macaca (Macaca mulatta), chimpanzee (Pan troglodytes)), treeshrew (Tupaia glis), mice (Mus musculus) and cow (Bos taurus). To identify positively selected genes of marmoset, the marmoset was set as foreground when detecting PSGs in the lineage, while the non-NWM were set as foreground when to identify NWM-specific PSGs. A total of 13,995 one-to-one orthologous genes were identified among ten species. To minimize the effect of gene annotation, we retrieved the corresponding CDS that shared the same isoform with human. These genes were used as an input dataset to conduct multiple sequence alignment using PRNAK (v.170427) and guidance (v2.02) to improve alignment. 





### Detection of positive selection genes

The positive selection sites within a specific lineage were detected by the branch-site model. Likelihood ratio test (LRT) was used to compare a model allowing sites to be under positive selection on the foreground branch with a null model in which sites may evolve neutrally and under purifying selection. The adjusted p values were also calculated using the FDR method, then genes with an adjusted p value less than 0.05 were treated as candidates for positive selection. The numbers of PSGs detected in marmoset and NWM, according to these criteria, were 324 and 59, respectively. To minimize effects of alignment, we filtered genes based on the condition of its positively selected sites following these criterions,
 
- sites with gap number more than 2 were excluded; 
- sites with nonsynonymous substitutions larger than 2 were excluded; 
- more complicated cases account for manual check. 

If one gene had no confident site, the gene would be removed. After filtering, the numbers of PSGs with high confidence detected in marmoset and NWM, were 207 and 38, respectively. The GO and KEGG pathway enrichment analysis were conducted using [Kobas3 website] (http://kobas.cbi.pku.edu.cn/kobas3) with Fisher's exact test and Benjamini and Hochberg method for multiple testing correction.
