# Hammock changelog
## Version 1.1.2
### 2017-11-29
- multi-line fasta files are now accepted as input
- new parameter: `--temp` to change the temporal files directory location from command line
- project now fully tracked on GitHub
- markdown documentation

## Version 1.1.1
### 2017-10-02
- parameter `-h, --min_match_states` has been renamed to `-h, --min_conserved_positions`, which better describes a slight change in its meaning. `--min_conserved_positions` defines the minimal number of each cluster's MSA positions that satisfy both `-k, --min_ic` and `-y, --max_gap_proportion`. The thing is, since version 1.1.0, HMM match states are defined as ALL the MSA positions between the leftmost and the rightmost conserved position (if inner gaps are not allowed). So in theory, with 2 conserved positions, there can be more than 2 match states. Since version 1.1.1, parameter `-h, --min_conserved_positions` defines the minimum nuber of conserved positions, not match states.
- New parameter `-as --additional_sequences` is available in "cluster" mode. It accepts a path to a fasta file containing sequences that will be added to the sequence pool during clustering. This is useful e.g. when additional rounds of cluster extension are required after a Hammock run has finished. In such a case, -as can be used to include `final_remaining_sequences.fa` file in a convenient way.  
- `-n, --assign_thresholds` and `-r, --merge_thresholds` sequences now may contain negative values. When a negative value is found, the corresponding extension/merging step is skipped. 
- Hammock now produces less temporal files. The difference is significant for very large datasets.
- In "cluster" mode, input files containing unaligned clusters are now allowed.
- Some low priority bugs fixed

## Version 1.1.0
### 2017-06-22
- New initial clustering algorithm (clinkage) for small datasets (up to 10 000 unique sequences by default) 
- New initial clustering algorithm (greedy clinkage) for large datasets (over 10 000 unique sequences by default). The original greedy clustering algorithm has been removed.
- Cluster cores now can be selected on the basis of their unique size (the `-U, --unique` switch)
- New options: `--min_correlation` `--min_cluster_size` `--min_cluster_unique_size`
- Options `-a, --part_threshold` and `-s, --size_threshold` are now considered deprecated due to confusion. Nevertheless, they will still work.
- New order option for greedy incremental clustering: `--order input`
- When `max_inner_gaps` is 0, HMM match states are defined as ALL the msa positions between the lefmost and the rightmost position satisfying `--min_ic` and `--max_gap_proportion`, even if some of the positions in between do not satisfy these two. Default score thresholds were slightly adjusted (increased) to reflect this
- When there are no database sequences left, Hammock continues to cluster by only performing cluster merging steps
- KLD is now calculated even when there are clusters only containing a single unique sequence (these are omitted form KLD calculation)
- Fixed a bug rarely causing incorrect alignments when sequences of variable lengts are present in the dataset
- From now on, we prefer questions and notes to be directed to GitHub: (<https://github.com/krejciadam/hammock/issues/>)

## Version1.0.6
### 2016-05-30
- Bugs causing problems when running in `cluster` mode fixed
- Hammock now produces (much) less temporal files when running on big datasets
- More logging and file existence controls

## Version 1.0.5
### 2016-03-18
- Fixed a bug causing problems when `-f tab` is in use
- Default values of `--max_shift` and `--min_match_states` are now guessed on the bases of average sequence length
- Updates in logging and error messages

## Version 1.0.4
### 2016-01-15
- all input/output tabular files are now tab-delimited
- `.cl` format no longer in use (`.cl` files are not saved and can't be loaded). Clusters can be loaded in `.tsv` format
- New output file saved for every run in greedy and full mode: `_clusters_sequences_original_order.tsv`
- Added controls for external tools: if they don't run properly, error messages should be more informative

## Version 1.0.3
### 2015-11-24
- Bug fixes
- KLD reporting should now work properly

## Version 1.0.2
### 2015-07-29
- Added new parameters: 
`-R, --order [size, alphabetic, random]`   The order of sequences during greedy clustering
`-S, --seed` A seed to make random processes deterministic (if `-R, --random` is in use)


## Version 1.0.1
### 2015-04-17
- Added the inline help option
- Bug fixes


## Version 1.0.0
- Initial version
