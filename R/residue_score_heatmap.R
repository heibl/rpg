#' Plot GUIDANCE Residue Scores and Alignment
#'
#' @param guidance_object is simply the output object of guidance
#' @param file if a dir path is supplied, the alignment is stored in a pdf, if file=NULL (default) alignment is plotted in R
#'
#' @author Franz-Sebastian Krah

guidance_heatmap <- function(guidance_object, file = NULL){
  txt <- as.vector(as.character(test$base_msa))
  GRSC <- data.frame(guidance_object$GUIDANCE_residue_score, txt)
  rown <- dim(guidance_obj$GUIDANCE_sequence_score)[1]
  coln <- dim(guidance_obj$GUIDANCE_residue_score)[1]/rown
  w <- coln/10
  h <- rown/4

  p <- ggplot(GRSC, aes(col, row)) +
    geom_tile(aes(fill = residue_score), colour = "white") +
    scale_fill_gradient(low = "red", high = "yellow")+
    scale_y_reverse()+
    ylab("Sequences (input order)")+
    xlab("Sites")+
    guides(fill=guide_legend(title="MSA\nconfidence\nscale"))+
    theme_bw()+
    theme(legend.position="left")+
    theme(plot.title = element_text(size = 20, face = "bold") ,
      legend.title=element_text(size=18) ,
      legend.text=element_text(size=14))
  p

  if(!is.null(file)){
    pdf(file, width = w, height = h)
    p <- p + geom_text(data=GRSC,aes(label = txt))
    print(p)
    dev.off()
  }

  if(is.null(file)){ p }
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
