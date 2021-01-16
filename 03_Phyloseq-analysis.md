Analyse des donées avec Phyloseq
================
Vincent Noah
14 janvier 2021

  - [Stratification des données en fonction de l’objet de
    l’étude.](#stratification-des-données-en-fonction-de-lobjet-de-létude.)
      - [Création de plusieurs objets.](#création-de-plusieurs-objets.)
  - [Création de l’objet ps.](#création-de-lobjet-ps.)
  - [Assignation des données avec dna grâce au package Biostrings pour
    manipuler les données
    biologiques.](#assignation-des-données-avec-dna-grâce-au-package-biostrings-pour-manipuler-les-données-biologiques.)
  - [Création de la table.](#création-de-la-table.)
  - [Visualisation de la prévalence](#visualisation-de-la-prévalence)
      - [Visualisation en Bar plot](#visualisation-en-bar-plot)

Tout d’abord nous allons charger les packages nécessaires et charger
l’image de l’analyse avec DADA2, afin de pouvoir continuer.

``` r
library(phyloseq)
library(ggplot2)
```

``` r
load("02_Data-analysis-with-DADA2_FinalEnv")
```

# Stratification des données en fonction de l’objet de l’étude.

## Création de plusieurs objets.

Avec Logan, on a passé des heures à stratifier correctement nos données,
mais nous n’avons pas réussi bien que nous soyons très proches de
réussir. On a donc laissé ce code pour pouvoir continuer une analyse,
bien que celui -ci ne soit pas exacte.

``` r
samples.out <- rownames(seqtab.nochim)
placard <- sapply(strsplit(samples.out, "p"), `[`, 2)
frigo <- (sapply(strsplit(samples.out, "f"), `[`, 3))
frigo <- substr(placard,1,10)
samdf <- data.frame(Placard=placard, Frigo=frigo)
samdf$Placard <- c("p1","p2","p3","p4","p5", "p6","p7","p8","p9","p10")
samdf$Frigo[samdf$Placard==10] <- c("f1","f2","f3","f4","f5", "f6","f7","f8","f9","f10")
rownames(samdf) <- samples.out
```

# Création de l’objet ps.

``` r
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
```

# Assignation des données avec dna grâce au package Biostrings pour manipuler les données biologiques.

``` r
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 2887 taxa and 20 samples ]
    ## sample_data() Sample Data:       [ 20 samples by 2 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 2887 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 2887 reference sequences ]

# Création de la table.

On observe qu’il y a de nombreuses protéobactéries, actinobactérie, et
firmicute. Cela semble logique car ces bactéries sont très abondantes,
et les objets de l’étude se trouvent sur des surfaces que l’on touche
souvent.

``` r
table(tax_table(ps)[, "Phylum"], exclude = NULL)
```

    ## 
    ##  Abditibacteriota   Acidobacteriota  Actinobacteriota    Armatimonadota 
    ##                 9                68               402                 7 
    ##      Bacteroidota  Bdellovibrionota  Campilobacterota       Chloroflexi 
    ##               298                25                 4                85 
    ##     Crenarchaeota     Cyanobacteria  Deferribacterota      Deinococcota 
    ##                 8               100                 1                26 
    ##      Dependentiae  Desulfobacterota   Elusimicrobiota Entotheonellaeota 
    ##                 1                 2                 1                 1 
    ##     Euryarchaeota    Fibrobacterota        Firmicutes    Fusobacteriota 
    ##                 1                 3               347                 8 
    ##   Gemmatimonadota Methylomirabilota       Myxococcota      Nitrospirota 
    ##                40                 2                39                 1 
    ##   Patescibacteria   Planctomycetota    Proteobacteria           RCP2-54 
    ##                 6                48               894                 3 
    ##     Spirochaetota       Sumerlaeota  Thermoplasmatota Verrucomicrobiota 
    ##                 1                 2                 1                37 
    ##              <NA> 
    ##               416

``` r
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))


plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

    ##               Phylum        1    2
    ## 1   Abditibacteriota 1.333333   12
    ## 2    Acidobacteriota 1.352941   92
    ## 3   Actinobacteriota 2.460199  989
    ## 4     Armatimonadota 1.000000    7
    ## 5       Bacteroidota 1.791946  534
    ## 6   Bdellovibrionota 1.480000   37
    ## 7   Campilobacterota 1.000000    4
    ## 8        Chloroflexi 1.329412  113
    ## 9      Crenarchaeota 1.625000   13
    ## 10     Cyanobacteria 2.610000  261
    ## 11  Deferribacterota 1.000000    1
    ## 12      Deinococcota 2.538462   66
    ## 13      Dependentiae 1.000000    1
    ## 14  Desulfobacterota 1.000000    2
    ## 15   Elusimicrobiota 1.000000    1
    ## 16 Entotheonellaeota 1.000000    1
    ## 17     Euryarchaeota 1.000000    1
    ## 18    Fibrobacterota 1.000000    3
    ## 19        Firmicutes 2.308357  801
    ## 20    Fusobacteriota 2.000000   16
    ## 21   Gemmatimonadota 1.150000   46
    ## 22 Methylomirabilota 1.000000    2
    ## 23       Myxococcota 1.102564   43
    ## 24      Nitrospirota 2.000000    2
    ## 25   Patescibacteria 1.000000    6
    ## 26   Planctomycetota 1.041667   50
    ## 27    Proteobacteria 2.437360 2179
    ## 28           RCP2-54 1.333333    4
    ## 29     Spirochaetota 1.000000    1
    ## 30       Sumerlaeota 1.000000    2
    ## 31  Thermoplasmatota 1.000000    1
    ## 32 Verrucomicrobiota 1.189189   44
    ## 33              <NA> 1.163462  484

# Visualisation de la prévalence

Ce graphique nous permet de confirmer cette prévalence. On observe que
ce sont seulement certain phyla qui sont très présent.

``` r
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

![](03_Phyloseq-analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Visualisation en Bar plot

``` r
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Frigo", fill="Family") + facet_wrap(~Placard, scales="free_x")
```

![](03_Phyloseq-analysis_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Frigo", fill="Family")
```

![](03_Phyloseq-analysis_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Placard", fill="Family") 
```

![](03_Phyloseq-analysis_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Quand on observe tous ces graphiques, on remarque que l’abondance des
bactéries est variable d’un placard à l’autre, d’un frigo à l’autre,
cependant, on peut observer certaines similarités. C’est par exemple ce
que l’on observe avec les Streptocoques. Tout comme logan, il y a la
présence de staphylocoques. Les staphylocoques tout comme les
Streptocoques, peuvent être pathogènes, bien que l’on retrouve ces
bactéries naturellement sur la peau en grande quantité. Contrairement à
Logan, il n’y a pas la présence d’entérobactéries. Ces bactéries sont
des marqueurs de contamination. Cette étude peut donc permettre de faire
prendre conscience en bien comme en mal de l’omniprésence des bactéries
sur de nombreuses surfaces que l’on utilise régulièrement.

Malheureusement, comme nous n’avions pas pu séparer correctement nos
données, on n’a pas pu réaliser une analyse très détaillée de la
présence des bactéries dans les cuisines. On aurait pu observer la
différence que l’on retrouve entre les frigos et les placards, mais
aussi, s’il y a de grandes différences entre les cuisines.
