# RDP-classifier-training-sets

The [RDP native training sets](http://rdp.cme.msu.edu/classifier/classifier.jsp;jsessionid=49C6531478DD4D70648EC6FC695B8EA3.10.0.0.9) were commonly employed for amplicon classification of bacteria and archaea. This was supplemented here with the Silva reference sequences, which can be used as a complement to the absence of taxonomic information if necessary. The latest Silva reference sequence sets [(v138)](https://www.arb-silva.de/no_cache/download/archive/release_138/), including bacteria, archaea and eukaryotes, have been reformatted for using the RDP classifier [(v2.13)](https://sourceforge.net/projects/rdp-classifier/). The latest(29.11.2022) [Unite + INSD](https://unite.ut.ee/repository.php) reference sequences were also trained to unify the workflow interface.

Training sets and reference sequences are available at [https://github.com/KaiMa-endeavour/RDP-classifier-training-sets/releases](https://github.com/KaiMa-endeavour/RDP-classifier-training-sets/releases).

## Training sets overview
|Sources|Taxonomic objects|Number of unique sequences|Ranks|
|:----:|:----:|:----:|:----:|
|native RDP|Bacteria|2833748|domain, phylum, class, order, family, genus|
|native RDP|Archaea|145117|domain, phylum, class, order, family, genus|
|Silva v138|Bacteria|1599148|domain, phylum, class, order, family, genus|
|Silva v138|Archaea|62300|domain, phylum, class, order, family, genus|
|Silva v138|Eukaryotes|149373|domain, kingdom, phylum, class, order, family, genus|
|Unite + INSD|Fungi|5335995|kingdom, phylum, class, order, family, genus, species|

## How to use
```
conda install -c bioconda rdp_classifier
rdp_classifier -Xmx16g classify -t ./training set/rRNAClassifier.properties -o ./rdp_output.txt ./query.fasta
```

The missing fields in the headers were filled with `NA` to align the lineages of each individual and to be compatible with the 'rdp_classifier train' function.

## References

Cole, J. R., Wang, Q., Fish, J. A., Chai, B., McGarrell, D. M., Sun, Y., Brown, C. T., Porras-Alfaro, A., Kuske, C. R., & Tiedje, J. M. (2014). Ribosomal Database Project: Data and tools for high throughput rRNA analysis. Nucleic Acids Research, 42(D1), D633–D642. https://doi.org/10.1093/nar/gkt1244

Nilsson, R. H., Larsson, K.-H., Taylor, A. F. S., Bengtsson-Palme, J., Jeppesen, T. S., Schigel, D., Kennedy, P., Picard, K., Glöckner, F. O., Tedersoo, L., Saar, I., Kõljalg, U., & Abarenkov, K. (2019). The UNITE database for molecular identification of fungi: Handling dark taxa and parallel taxonomic classifications. Nucleic Acids Research, 47(D1), D259–D264. https://doi.org/10.1093/nar/gky1022

Quast, C., Pruesse, E., Yilmaz, P., Gerken, J., Schweer, T., Yarza, P., Peplies, J., & Glöckner, F. O. (2013). The SILVA ribosomal RNA gene database project: Improved data processing and web-based tools. Nucleic Acids Research, 41(D1), D590–D596. https://doi.org/10.1093/nar/gks1219

Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Applied and Environmental Microbiology, 73(16), 5261–5267. https://doi.org/10.1128/AEM.00062-07
