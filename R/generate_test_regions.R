


#' Generate set 1 of testing regions.
#'
#' @return Nothing, saves files in the inst/testdata folder
#' @export
#'
#' @examples
#'
#' # no run
#'
generate_test_regions <- function() {

  set.seed(42)

  ofolder = "inst/testdata"
  fs::dir_create(ofolder)

  library(GenomicRanges)
  library(magrittr)

  g = helperMut::genome_selector()
  rg = helperMut::get_region_chunks(gr = g,wl = 1000000,unlist = TRUE)
  rg2 = helperMut::get_region_chunks(gr = g,wl = 505002,unlist = TRUE)
  rg3 = helperMut::get_region_chunks(gr = g,wl = 696165 ,unlist = TRUE)

  rg = rg[seqnames(rg) %in% c("chr1","chr2")]
  rg2 = rg2[seqnames(rg2) %in% c("chr1","chr2")]
  rg3 = rg3[seqnames(rg3) %in% c("chr1","chr2")]

  mcols(rg)$bin = sample(x = letters[1:3],size = length(rg),replace = TRUE)
  mcols(rg2)$bin = sample(x = LETTERS[11:13],size = length(rg2),replace = TRUE)
  mcols(rg3)$bin = sample(x = "S",size = length(rg3),replace = TRUE)

  rg_list = GenomicRanges::split(rg,f = mcols(rg)$bin)
  rg2_list = GenomicRanges::split(rg2,f = mcols(rg2)$bin)
  rg3_list = GenomicRanges::split(rg3,f = mcols(rg3)$bin)

  lapply(rg_list, function(x){
    #browser()
    reg_reduced = reduce(x)
    unique_bin = mcols(x)$bin %>% unique
    oname = fs::path(ofolder,glue::glue("test_bed1_bin{unique_bin}.bed"))
    rtracklayer::export.bed(object = reg_reduced,con = oname)
  })

  lapply(rg2_list, function(x){
    reg_reduced = reduce(x)
    unique_bin = mcols(x)$bin %>% unique
    oname = fs::path(ofolder,glue::glue("test_bed2_bin{unique_bin}.bed"))
    rtracklayer::export.bed(object = reg_reduced,con = oname)
  })

  # this will have strand bias
  lapply(rg3_list, function(x){


    unique_bin = mcols(x)$bin %>% unique

    strand_tmp = sample(c("+","-"),
                        size =  length(x),
                        replace = T)

    oname = fs::path(ofolder,glue::glue("test_bed3_bin{unique_bin}_featOnPlus.bed"))
    out_obj = x[strand_tmp == "+"]
    reg_reduced = reduce(out_obj)
    rtracklayer::export.bed(object = reg_reduced,con = oname)

    oname = fs::path(ofolder,glue::glue("test_bed3_bin{unique_bin}_featOnMinus.bed"))
    out_obj = x[strand_tmp == "-"]
    reg_reduced = reduce(out_obj)
    rtracklayer::export.bed(object = reg_reduced,con = oname)
  })

  master_name = c(rep("bed1",3),rep("bed2",3))
  ref_letters = c(letters[1:3],LETTERS[11:13])

  ref_table = data.frame(
    master_name = master_name,
    bins = ref_letters,
    file_name = fs::path(ofolder,glue::glue("test_{master_name}_bin{ref_letters}.bed")),
    strand = "*"
  )

  sfname = c(
    "test_bed3_binS_featOnPlus.bed",
    "test_bed3_binS_featOnMinus.bed"
  )

  ref_table_s = data.frame(
    master_name = c("bed3","bed3"),
    bins = c("S_plus","S_minus"),
    file_name = fs::path(ofolder,sfname),
    strand = c("+","-")
  )
  readr::write_tsv(x = rbind(ref_table),
                   path = "inst/testdata/test_bins_nostrand.tsv")
  readr::write_tsv(x = rbind(ref_table,ref_table_s),
                   path = "inst/testdata/test_bins.tsv")


  # Now I generate the mutations file
  library(VariantAnnotation)

  nmuts = c(1500,1000)
  total_size = seqlengths(g)[1:2]

  pos1 = sample(size = nmuts[1],x = total_size[1],replace = F)
  pos2 = sample(size = nmuts[2],x = total_size[2],replace = F)

  VRanges(
    seqnames = c(rep("chr1", nmuts[1]),rep("chr2", nmuts[2])),
    ranges = IRanges(
      start = c(pos1,pos2),
      width = 1
    ),
    ref = "R",
    alt = "A",
    sampleNames = "sampleA"
  ) -> dat_vr
  dat_vcf = as(object = dat_vr,Class = "VCF")
  writeVcf(obj = dat_vcf,filename = fs::path(ofolder,"test_vcf.vcf"))

}


#' Generate set 2 of testing regions.
#'
#' @return Nothing, saves files in the inst/testdata folder
#' @export
#'
#' @examples
#'
#' # no run
#'
generate_test_regions2 <- function() {
  set.seed(42)

  library(GenomicRanges)
  library(magrittr)


  path = "inst/testdata/"

  channel1 = GRanges(
    seqnames = "chr1",
    ranges = IRanges(
      start = c(1,100,200) +  1e6,
      width = 100
    ),
    id = c("a","b","c")
  )

  channel1 %>% base::split(.$id) %>%
    purrr::walk(function(x){
      uid = unique(x$id)
      rtracklayer::export.bed(object = x,
                              con = fs::path(path,
                            glue::glue("channel1_bin{uid}.bed")))
    })

  channel2 = GRanges(
    seqnames = "chr1",
    ranges = IRanges(
      start = c(1,200,50,150,250,100) +  1e6,
      width = 50
    ),
    id = c("K","K","L","L","L","M")
  )

  channel2 %>% base::split(.$id) %>%
    purrr::walk(function(x){
      uid = unique(x$id)
      rtracklayer::export.bed(object = x,
                              con = fs::path(path,
                                   glue::glue("channel2_bin{uid}.bed")))
    })

  channel3 = GRanges(
    seqnames = "chr1",
    ranges = IRanges(
      start = c(1,150) +  1e6,
      width = 150
    ),
    id = c("S","S")
  )

  channel3 %>% base::split(1:2) %>%
    purrr::walk2(c("plus","minus"),function(x,name){
      uid = unique(x$id)
      rtracklayer::export.bed(object = x,
                              con = fs::path(path,
                     glue::glue("channel3_{uid}_{name}.bed")))
    })

  table_tsv = data.frame(
    master_name = c(rep("Channel1",3),rep("Channel2",3),rep("Channel3",2)),
    bins = c(c("bina","binb","binc"),c("binK","binL","binM"),"S_plus","S_minus"),
    strand = c(rep("*",6), "+","-")
  )

  table_tsv$file_name = fs::path(path,
          glue::glue("{table_tsv$master_name}_{table_tsv$bins}.bed"))

  readr::write_tsv(x = table_tsv,path = fs::path(path,"channels_bins.tsv"))


  library(VariantAnnotation)

  dat_vr = VRanges(
    seqnames = "chr1",
    ranges = IRanges(
      start = c(25,75,125,175,225,275) + 1e6,
      width = 1
    ),
    ref = "N", # they are used to check for indels only I think
    alt = "H",
    sampleNames = "David"
  )

  as(object = dat_vr,Class = "VCF") -> dat_vcf

  writeVcf(dat_vcf,filename = fs::path(path,"david_muts.vcf"))

}

