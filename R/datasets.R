
#' Conversion of the CIS-BP motif database to k-mer space.
#'
#' A list containing: logical matrices of cognate (T) and non-cognate (F) k-mers for each human motif in CIS-BP (as defined in de Boer and Regev, 2017, bioRxiv);
#' information on how these matrices were created; information on the motifs themselves; information on the similarity between motifs.
#'
#' $binaryPWMScores a matrix of k-mers by motifs, where entry is True if k-mer is cognate of motif, and False otherwise; made from CIS-BP PWM motifs
#' $binaryPBMZScores a matrix of k-mers by motifs, where entry is True if k-mer is cognate of motif, and False otherwise; made from CIS-BP PBM Z-scores
#' $info a list of information describing how this was created.
#' $TFTable a data.frame with information about each motif (from CIS-BP: \url{http://cisbp.ccbr.utoronto.ca/})
#' $similarMotifs a data.frame containing Motif_IDs and the Pearson r between the two (for those with r>0.5)
#' @source \url{http://cisbp.ccbr.utoronto.ca/}
"cisbp"
