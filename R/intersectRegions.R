


#' Read input regions
#'
#' It needs a design table with the files used as regions of interest.
#' It then reads and structures the bed files accordingly.
#'
#' Used in regionMut region subcommand.
#'
#' @param ref_data a dataframe with the design table
#'
#' @return A list of list with the regions of interest.
#' @export
#'
#' @examples
#'
#' # no run
#'
read_input_regions <- function(ref_data) {

  ref_data_strand = ref_data[ref_data$strand != "*",]
  #ref_data = ref_data[ref_data$strand == "*",]

  if (nrow(ref_data_strand)> 2){
    # I think if there are more than 2 strand channels then the problem of
    # one refering to one of the strands and the other to the reverse is
    # problematic. Basically the ntps at risk, and the offset term will be
    # wrong as we would need to count the same bp twice? (I am still not sure)
    stop("More than 2 strand channels not supported")
  }

  # reading regions ========
  # this reads all the files
  lol_input = ref_data %>% base::split(.$master_name) %>%
    purrr::map(function(x){
      GRs = x$file_name %>% as.character() %>%
        purrr::map(rtracklayer::import.bed) %>%
        purrr::map2(x$strand,
                    function(x,st){
                      strand(x) = st
                      x
                    })
      tm_l = GRs
      names(tm_l) = x$bins
      return(tm_l)
    })

  return(lol_input)

}

#' Intersect regions
#'
#' It performs a combinatorial intersection of regions. Basically, it
#' combines features in the different channels by intersectling them
#' with all the available combinations
#'
#' @param list_of_regions The list of regions obtained from \code{read_input_regions}
#'
#' @return a list of intersected regions
#' @export
#'
#' @examples
#'
#' # no run
process_region_interactions2 <- function(list_of_regions){

    # this are all the bins in each region channel
    bins_per_feature = purrr::map(list_of_regions,names)
    comb = expand.grid(bins_per_feature)

    # make intersections ==========
    id = 1:nrow(comb)

    # this is slow
    id %>% purrr::map(function(i){

      GRs = list()
      for (j in 1:ncol(comb)){

        master_name = colnames(comb)[j]

        GRs[j] = list_of_regions[[master_name]][[comb[i,master_name]]]
      }

      GRs %>% purrr::map(strand) %>% purrr::map(unique) %>% unlist -> stds
      unique_strand = stds[stds != "*"]
      stopifnot(length(unique_strand) <= 1)

      if (length(unique_strand) == 0){
        # this means that there is no strand in the input file
        intersect_regions = gr_intesect(gr_base_list = GRs)
      } else {
        # when at least there is some strands
        intersect_regions = gr_intesect(gr_base_list = GRs,
                                        ignore.strand = TRUE)
        strand(intersect_regions) = unique_strand
      }

      if(length(intersect_regions) == 0){
        mssg_war = glue::glue(
          "Interaction was not possible at id: {i}"
        )
        warning(mssg_war)
        return(GRanges())
      } else{
        MCOLS = mcols(intersect_regions)
        MCOLS = cbind(MCOLS,comb[i,])
        MCOLS$id = i
        mcols(intersect_regions) = MCOLS
        return(intersect_regions)
      }
    }) -> regions_list

    return(regions_list)

  }



#' Oligo Counts
#'
#' Computes the olgonucleotide counts from GRanges object.
#'
#' @param gr a GRanges object
#' @param genome a genome object from class BSgenome
#' @param k Number of sites to expand the oligonucleotide, size = 2*k + 1
#' @param refSet A reference set of positions.
#' @param data.frame boolean if to return a dataframe or a vector
#'
#' @return
#'
#' A dataframe (or a vector) with the oligonucleotides contined in the
#' selected region.
#'
#' @export
#'
#' @examples
#'
#' library(GenomicRanges)
#' gr = GRanges(seqnames = c("chr1"),
#'         ranges = IRanges(start = c(1e6,2e6),width = 100))
#' genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#'
#' obtain_oligo_counts(gr = gr,
#'                     genome = genome,
#'                     k = 1,
#'                     refSet = c("C","A"),
#'                     data.frame = TRUE)
#'
obtain_oligo_counts <- function(gr,
                                genome,
                                k,
                                refSet,
                                data.frame = FALSE) {
  size = (2*k)+1
  # so the function getSeq needs us to provide also the
  # extra basepair (or k) that it's missing from the
  # ranged input
  gr = helperMut::extend(x = gr,
                         upstream = k,
                         downstream = k)

  genome_sql = seqlengths(genome)
  seqlengths(gr) = suppressWarnings(genome_sql[seqlevels(gr)])
  gr = GenomicRanges::trim(gr)
  reg = Biostrings::getSeq(genome,gr)
  oligo_counts = Biostrings::oligonucleotideFrequency(x = reg,
                                                      width = size)
  oligo_counts = apply(oligo_counts,2,sum)

  nnames = helperMut::simplify_ctx(ctx = names(oligo_counts),
                                   simplify_set = refSet)
  unlist(
    lapply(
      base::split(oligo_counts,
                  nnames),
      sum)
  ) -> oligo_counts_vec

  if(data.frame){
    data.frame(
      ctx = names(oligo_counts_vec),
      ctx_counts = oligo_counts_vec
    )
  } else {
    oligo_counts_vec
  }
}

#' Get offset
#'
#' Obtain the nucleotides at risk of each region.
#'
#' It inputs the regions from \code{process_region_interactions2}.
#'
#' @param gr_list a list of GRanges objects
#' @param k Number of sites to expand the oligonucleotide, size = 2*k + 1
#' @param genome a genome object from class BSgenome
#' @param refSet A reference set of positions.
#'
#' @return A data frame with oligonucleotide counts to use as offset.
#' @export
#'
#' @examples
#'
#' # not run
#'
get_offset <- function(gr_list,
                       k,
                       genome,
                       refSet
) {

  ## here I am introducng the factors
  if (is.null(names(gr_list))){
    names(gr_list) = factor(as.character(1:length(gr_list)))
  }

  lapply(gr_list,
         obtain_oligo_counts,
         genome = genome,
         k = k,
         refSet = refSet,
         data.frame = TRUE) -> ocounts_df_list

  res_df = dplyr::bind_rows(ocounts_df_list,.id = "combination_id")

  return(res_df)
}



# utils -------------------------------------------------------------------

#' Intersect a list of GRanges
#'
#' Given a list of GRanges objects it returns the intersection of all of them.
#' Thus, the resulting sites represent sites that are shared in all the
#' GRanges instances in the list.
#'
#' @param gr_base_list list with GRanges
#' @param ignore.strand If to ignore strand when intersecting regions.
#'
#' @return a single GR object with the intersected sites
#' @export
#'
#' @examples
#'
#' # no run.
#'
gr_intesect <- function(gr_base_list,ignore.strand = FALSE){

  base_gr = gr_base_list[[1]]

  # this handles the case when only on GR object is inputed
  if (length(gr_base_list) > 1){
    # this keeps adding the rest of the GRs in the base object
    for (i in 2:length(gr_base_list)){

      current_gr = gr_base_list[[i]]

      base_gr = GenomicRanges::intersect(base_gr,
                                         current_gr,
                                         ignore.strand = ignore.strand)

    }
  }
  return(base_gr)
}
