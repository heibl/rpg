#' Plot GUIDANCE Residue Scores
#' @description Plots a heatmap of GUIDANCE residue scores
#' @param obj is the output object of guidance
#' @param file if a path is supplied, the alignment is stored in a pdf, if file=NULL (default) alignment is promted in R
#'
#' @author Franz-Sebastian Krah
#' @import ggplot2
#' @import scales
#' @export

heatmap.msa <- function(obj, file = NULL){

  txt <- as.vector(as.character(obj$base_msa))
  if (length(grep("scores", names(obj))) > 0)
    obj <- obj$scores

  mat <- data.frame(obj$residue_pair_residue_score, txt)
  rown <- max(mat$residue)
  coln <- max(mat$col)
  res_mat <- matrix(mat$score, nrow = rown, ncol = coln)

  ## scales for PDF
  w <- coln/10
  h <- rown/4

  p <- ggplot(mat,aes(col, residue)) +
    geom_tile(aes(fill = score), colour = "white") +
    scale_fill_gradient(low = "red", high = "yellow") +
    ylab("Sequences (input order)") +
    xlab("Sites") +
    scale_y_reverse(breaks= pretty_breaks()) +
    theme(plot.title = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 14))+
    ylab("Sequences (input order)")+
    xlab("Sites")+
    theme_bw()+
    theme(legend.position="left")+
    theme(legend.title=element_text(size=18) ,
      legend.text=element_text(size=12),
      legend.position="top",
      legend.direction="horizontal")+
    guides(fill=guide_legend(title="MSA confidence scale (1=high)"))

  if (!is.null(file)){
    pdf(file, width = w, height = h)
    p <- p + geom_text(data = mat, aes(label = txt))
    print(p)
    dev.off()
  } else {
    return(p)
  }
}
