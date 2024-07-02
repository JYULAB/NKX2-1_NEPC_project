library(epiAneufinder)
epiAneufinder(input="/projects/p20023/Viriya/analysis/foxa2/multiome/arc_ranger/d2/atac_fragments.tsv.gz", 
              outdir="./", 
              blacklist="../../hg38-blacklist.v2.bed", 
              windowSize=1e6, 
              genome="BSgenome.Hsapiens.UCSC.hg38", 
              exclude=c('chrX','chrY','chrM'), 
              reuse.existing=TRUE,
              title_karyo="Karyogram of D2 data", 
              ncores=6,
              minFrags=10000,
              minsizeCNV=0,
              k=4,
              plotKaryo=TRUE)

epiAneufinder(input="/projects/p20023/Viriya/analysis/foxa2/multiome/arc_ranger/d14/atac_fragments.tsv.gz", 
              outdir="./", 
              blacklist="../../hg38-blacklist.v2.bed", 
              windowSize=1e6, 
              genome="BSgenome.Hsapiens.UCSC.hg38", 
              exclude=c('chrX','chrY','chrM'), 
              reuse.existing=TRUE,
              title_karyo="Karyogram of D14 data", 
              ncores=6,
              minFrags=10000,
              minsizeCNV=0,
              k=4,
              plotKaryo=TRUE)

epiAneufinder(input="/projects/p20023/Viriya/analysis/foxa2/multiome/arc_ranger/d21_add/atac_fragments.tsv.gz", 
              outdir="./", 
              blacklist="../../hg38-blacklist.v2.bed", 
              windowSize=1e6, 
              genome="BSgenome.Hsapiens.UCSC.hg38", 
              exclude=c('chrX','chrY','chrM'), 
              reuse.existing=TRUE,
              title_karyo="Karyogram of D21 Add data", 
              ncores=6,
              minFrags=10000,
              minsizeCNV=0,
              k=4,
              plotKaryo=TRUE)




