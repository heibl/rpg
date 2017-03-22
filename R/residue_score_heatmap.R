#' @title Plot GUIDANCE Residue Scores
#' @description Plots a heatmap of GUIDANCE residue scores
#' @param obj is the output object of guidance
#' @param file if a path is supplied, the alignment is stored in a pdf, if file=NULL (default) alignment is promted in R
#'
#' @author Franz-Sebastian Krah
#' @importFrom ggplot2 ggplot
#' @export

heatmap.msa <- function(obj, file = NULL){

  hot <- grep("reliability", names(obj)[1])
  guid <- grep("confidence", names(obj)[1])

  txt <- as.vector(as.character(obj$base_msa))

  if(length(guid) > 1){
    mat <- data.frame(obj$GUIDANCE_residue_score, txt)
    rown <- dim(obj$GUIDANCE_sequence_score)[1]
    coln <- dim(obj$GUIDANCE_residue_score)[1]/rown
  }
  if(length(hot) > 1){
    mat <- data.frame(obj$residue_reliability, txt)
    rown <- dim(obj$sequence_reliability)[1]
    coln <- dim(obj$column_reliability)[1]/rown
  }

  w <- coln/10
  h <- rown/4

  p <- ggplot(mat,aes(col, row)) +
    geom_tile(aes(fill = residue_score), colour = "white") +
    scale_fill_gradient(low = "red", high = "yellow") +
    ylab("Sequences (input order)") +
    xlab("Sites") +
    theme(plot.title = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 14))+
    scale_y_reverse()+
    ylab("Sequences (input order)")+
    xlab("Sites")+
    theme_bw()+
    theme(legend.position="left")+
    theme(legend.title=element_text(size=18) ,
      legend.text=element_text(size=14),
      legend.position="top",
      legend.direction="horizontal")+
    guides(fill=guide_legend(title="MSA confidence scale (1=high)"))

  if (!is.null(file)){
    pdf(file, width = w, height = h)
    p <- p + geom_text(data = GRSC, aes(label = txt))
    print(p)
    dev.off()
  } else {
    return(p)
  }
}



# guidance_score_heatmap <- function(guidance_obj, gap.col ="white"){
#   rown <- dim(guidance_obj$GUIDANCE_sequence_score)[1]
#   coln <- dim(guidance_obj$GUIDANCE_residue_score)[1]/rown
#   msa <- matrix(guidance_obj$GUIDANCE_residue_score, nrow = rown, ncol = coln)
#   hmcol<-brewer.pal(11,"RdYlBu")
#   require("gplots")
#   heatmap.2(msa,
#     # if(print.base)
#     # cellnote=mat2,
#     dendrogram = "none",
#     labRow = labels(guidance_obj$base_msa),
#     labCol = "",
#     Rowv = FALSE,
#     Colv = FALSE,
#     notecex=1.0,
#     notecol="black",
#     na.color=gap.col,
#     trace = "none",
#     key = T,
#     key.title = "MSA confidence\nscale",
#     key.ylab = "",
#     key.xlab = "Uncertain <---> Conficent",
#     density.info = "density",
#     keysize = 1.7,
#     col  = hmcol)
# }
#   if(is.null(file)){
#     print(p)
#   }
# }
