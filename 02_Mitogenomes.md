02_mitogenomes
================
2025-11-12

We aligned ALL the files to the mitochondrial genome which is from the
same animal as the nuclear genome we used as a reference.

Ran scripts 00.stagein.sh to get all files from Acacia. Then
01.fix_bam_headers.sh, 02.sort_bams.sh and 03.call_hapl_consensus.sh.
Next up, we converted to fasta using 04.convert_to_fasta.sh. Then I
downloaded goalign and ran 05.goalign.sh

I used an online tool to convert fasta to nexus format.

Then I ran 07.run_iqtree.sh and downloaded the resulting .contree file
in here.

``` r
tree <- read.tree("02_Mitogenomes_local//AllHapsClean_all.fasta.contree")
# Root the tree using sample "628131"
rooted_tree <- root(tree, outgroup = "628131", resolve.root = TRUE)

# Plot it
plot(rooted_tree, show.tip.label = TRUE)
```

![](02_Mitogenomes_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# Save to file
write.tree(rooted_tree, file = "02_Mitogenomes_local//AllHapsClean_rooted.tree")
```

    ## <ggproto object: Class ScaleDiscrete, Scale, gg>
    ##     aesthetics: colour
    ##     axis_order: function
    ##     break_info: function
    ##     break_positions: function
    ##     breaks: Shark Bay Exmouth Gulf Barrow Island Port Hedland Broome ...
    ##     call: call
    ##     clone: function
    ##     dimension: function
    ##     drop: TRUE
    ##     expand: waiver
    ##     get_breaks: function
    ##     get_breaks_minor: function
    ##     get_labels: function
    ##     get_limits: function
    ##     get_transformation: function
    ##     guide: legend
    ##     is_discrete: function
    ##     is_empty: function
    ##     labels: waiver
    ##     limits: function
    ##     make_sec_title: function
    ##     make_title: function
    ##     map: function
    ##     map_df: function
    ##     n.breaks.cache: NULL
    ##     na.translate: TRUE
    ##     na.value: grey50
    ##     name: waiver
    ##     palette: function
    ##     palette.cache: NULL
    ##     position: left
    ##     range: environment
    ##     rescale: function
    ##     reset: function
    ##     train: function
    ##     train_df: function
    ##     transform: function
    ##     transform_df: function
    ##     super:  <ggproto object: Class ScaleDiscrete, Scale, gg>

![](02_Mitogenomes_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](02_Mitogenomes_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

Focusing only on Australia now:
![](02_Mitogenomes_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

**NO CONTROL REGION**

I first downloaded seqkit: wget
<https://github.com/shenwei356/seqkit/releases/download/v2.3.1/seqkit_linux_amd64.tar.gz>

Then I had to index the reference

Then I had to adjust a lot of things, check scripts 08 - 11

![](02_Mitogenomes_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

![](02_Mitogenomes_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Now, I’m gonna try to run a tree with the full sequences. I am using the
online mafft tool to make an alignment.

![](02_Mitogenomes_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->![](02_Mitogenomes_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->
