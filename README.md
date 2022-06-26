# RDP-classifier-training-sets

The [RDP native training sets](http://rdp.cme.msu.edu/classifier/classifier.jsp;jsessionid=49C6531478DD4D70648EC6FC695B8EA3.10.0.0.9) are typically used for bacteria and archaea, however, the latest Silva reference sequences can be used as a supplement in cases where taxonomic information is not ideal. The latest Silva reference sequence sets [(v138)](https://www.arb-silva.de/no_cache/download/archive/release_138/), including bacteria, archaea and eukaryotes, has been reformatted for using the RDP classifier [(v2.13)](https://sourceforge.net/projects/rdp-classifier/). The latest version of the [Unite](https://github.com/terrimporter/UNITE_ITSClassifier) database was used for the fungi ITS training set.

Training sets and reference sequences are available at [https://github.com/KaiMa-endeavour/RDP-classifier-training-sets/releases](https://github.com/KaiMa-endeavour/RDP-classifier-training-sets/releases). 

## Training sets overview
|source|domain|number of unique sequences|ranks|
|:----:|:----:|:----:|:----:|
|native RDP|Bacteria|2833748|domain, phylum, class, order, family, genus|
|native RDP|Archaea|145117|domain, phylum, class, order, family, genus|
|Silva v138|Bacteria|1599148|domain, phylum, class, order, family, genus|
|Silva v138|Archaea|62300|domain, phylum, class, order, family, genus|
|Silva v138|Eukaryotes|149373|domain, kingdom, phylum, class, order, family, genus|
|Unite|Fungi|733093|domain, phylum, class, order, family, genus|

## How to use
```
conda install -c bioconda rdp_classifier
rdp_classifier -Xmx16g classify -f fixrank -t ./training set/rRNAClassifier.properties -o ./rdp_output.txt ./query.fasta
```
In order to meet the requirements of the function `rdp_classifier train`, the missing fields are filled with `NA`, so it is recommended to delete the `NA` in rdp_output.txt.

## Reference

Cole, J. R., Wang, Q., Fish, J. A., Chai, B., McGarrell, D. M., Sun, Y., Brown, C. T., Porras-Alfaro, A., Kuske, C. R., & Tiedje, J. M. (2014). Ribosomal Database Project: Data and tools for high throughput rRNA analysis. Nucleic Acids Research, 42(D1), D633–D642. https://doi.org/10.1093/nar/gkt1244

Nilsson, R. H., Larsson, K.-H., Taylor, A. F. S., Bengtsson-Palme, J., Jeppesen, T. S., Schigel, D., Kennedy, P., Picard, K., Glöckner, F. O., Tedersoo, L., Saar, I., Kõljalg, U., & Abarenkov, K. (2019). The UNITE database for molecular identification of fungi: Handling dark taxa and parallel taxonomic classifications. Nucleic Acids Research, 47(D1), D259–D264. https://doi.org/10.1093/nar/gky1022

Quast, C., Pruesse, E., Yilmaz, P., Gerken, J., Schweer, T., Yarza, P., Peplies, J., & Glöckner, F. O. (2013). The SILVA ribosomal RNA gene database project: Improved data processing and web-based tools. Nucleic Acids Research, 41(D1), D590–D596. https://doi.org/10.1093/nar/gks1219

Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Applied and Environmental Microbiology, 73(16), 5261–5267. https://doi.org/10.1128/AEM.00062-07

