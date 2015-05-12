#packages needed to convert gene symbols to ensembl ids
library(AnnotationDbi)
library(org.Hs.eg.db)

library(rvest) #web scraper by hadley wickham
library(httr) #http package to send GET/POST http methods
library(magrittr) #useful pipe operator

# SNiPA URIs --------------------------------------------------------------

`%+%` <- paste0 #poor man's homemade concat operator

#define URI's that will be used for HTTP requests
snipa.base              <-"http://snipa.helmholtz-muenchen.de/snipa"
snipa.pairwise.page.uri <- snipa.base %+% "/index.php?task=pairwise_ld"
snipa.ld.form.uri       <- snipa.base %+% "/backend/snipaProxySearch.php"
snipa.rand.id.uri       <- snipa.base %+% "/backend/snipaTempdir.php"
snipa.get.status.uri    <- function(i) snipa.base %+% "/tmpdata/" %+% i %+% "/status.txt"
snipa.get.report.uri    <- function(i) snipa.base %+% "/tmpdata/" %+% i %+% "/report.txt"
snipa.get.result.uri    <- function(i) snipa.base %+% "/tmpdata/" %+% i %+% "/proxySearch.results.csv"
user.agent <- "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_2) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/42.0.2311.50 Safari/537.36"


# Helper Functions --------------------------------------------------------

#generate a random ID through Snipa server method
#this is used throughout all snipa operations
snipa.get.rand <- function() {
  s <- html_session(snipa.rand.id.uri,
                    user_agent(user.agent))
  stop_for_status(s)
  s %>% html() %>% html_node('p') %>% html_text()
}

#generate a list which is converted to urlencoded POST() parameters
snipa.generate.ld.form.params <- function(snps_sentinels=NULL,
                                          id=NULL,
                                          genomerelease='grch37',
                                          referenceset='1kgpp3v5',
                                          population='eur',
                                          annotation='ensembl77',
                                          snps_input_type='snps',
                                          snps_gene=NULL,
                                          snps_region_chr=NULL,
                                          snps_region_begin=NULL,
                                          snps_region_end=NULL,
                                          rsquare=0.8,
                                          incl_sentinel=1,
                                          incl_funcann=0,
                                          download=1,
                                          dyn_tables=0,
                                          pairwise=1){
  params <- mget(names(formals()))
  params$id <- ifelse(is.null(id), snipa.get.rand(), id)
  params
}

#submits given query through HTTP POST method to the given URI
snipa.submit.query <- function(params, uri) {
  s <- html_session(snipa.pairwise.page.uri)
  stop_for_status(s)
  scook <- cookies(s)

  res <- POST(uri,
              body=params,
              encode='form',
              user_agent(user.agent),
              do.call(set_cookies, scook))

  stop_for_status(res)
  res <- httr::content(res, type='application/json')
  if(res$errmessage != '')
    stop(res$errmessage)

  stopifnot(res$stepnum == res$totalstepnum)

  res <- GET(snipa.get.result.uri(params$id),
             user_agent(user.agent),
             do.call(set_cookies, scook))
  stop_for_status(res)
  res <- httr::content(res, as='text')

  read.delim(textConnection(res), stringsAsFactors = F)
}

# Main Functions ----------------------------------------------------------

#get LD information for given SNPs
snipa.get.ld.by.snp <- function(snps, #vector of sentinel SNPs
                                rsquare=0.8, #rsquared threshold
                                pairwise=F,  #pairwise=T returns pairwise LD values for given sentinels
                                annotation=F, #whether functional annotation of SNPs will be returned
                                population=c('eur', 'afr', 'amr', 'eas', 'sas'),
                                ...){

  snps <- paste(snps, collapse = '\n')

  #translate arguments to real form variables
  params <- snipa.generate.ld.form.params(snps_sentinels = snps,
                                          rsquare=rsquare,
                                          pairwise=as.integer(pairwise),
                                          incl_funcann = as.integer(annotation),
                                          population=match.arg(population),
                                          ...)
  snipa.submit.query(params, snipa.ld.form.uri)
}

snipa.get.ld.by.region <- function(chr, begin, end,
                                   rsquare=0.8, #rsquared threshold
                                   annotation=F, #whether functional annotation of SNPs will be returned
                                   population=c('eur', 'afr', 'amr', 'eas', 'sas'),
                                   ...){
  #strip 'chr' prefix off, if exists
  if(substr(chr, 1, 3)=='chr') chr <- substr(chr, 4, nchar(chr))

  #translate arguments to real form variables
  params <- snipa.generate.ld.form.params(snps_region_chr = chr,
                                          snps_region_begin = begin,
                                          snps_region_end = end,
                                          snps_input_type='region',
                                          rsquare=rsquare,
                                          pairwise=0,
                                          incl_funcann = as.integer(annotation),
                                          population=match.arg(population),
                                          ...)
  snipa.submit.query(params, snipa.ld.form.uri)
}

#get LD information for given gene
snipa.get.ld.by.gene <- function(gene,
                                 rsquare=0.8,
                                 annotation=F, #whether functional annotation of SNPs will be returned
                                 gene.id=c('ENSEMBL', 'SYMBOL', 'ENTREZID'),
                                 population=c('eur', 'afr', 'amr', 'eas', 'sas'),
                                 ...){

  gene.id <- match.arg(gene.id)

  #convert symbol/entrezid to ensembl id
  convertIDs <- function(ids, from) {
    suppressWarnings(selRes <- AnnotationDbi::select(org.Hs.eg.db,
                                                     keys=ids,
                                                     keytype=from,
                                                     columns=c(from,'ENSEMBL')))
    return(selRes[match(ids, selRes[, 1] ), 2])
  }

  if(gene.id != 'ENSEMBL') gene <- convertIDs(gene, gene.id)

  #translate arguments to real form variables
  params <- snipa.generate.ld.form.params(snps_gene = gene,
                                          snps_input_type = 'gene',
                                          rsquare=rsquare,
                                          pairwise=0,
                                          incl_funcann = as.integer(annotation),
                                          population=match.arg(population),
                                          ...)
  snipa.submit.query(params, snipa.ld.form.uri)
}
