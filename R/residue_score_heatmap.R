#' Plot GUIDANCE Residue Scores and Alignment
#'
#' @param guidance_object is simply the output object of guidance
#' @param file if a dir path is supplied, the alignment is stored in a pdf, if file=NULL (default) alignment is plotted in R
#'
#' @author Franz-Sebastian Krah

guidance_heatmap <- function(guidance_obj, file = NULL){

  txt <- as.vector(as.character(guidance_obj$base_msa))
  GRSC <- data.frame(guidance_obj$GUIDANCE_residue_score, txt)
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
    guides(fill=guide_legend(title="MSA confidence scale\n 1=high confidence")) +
    theme_bw()+
    theme(legend.title=element_text(size=18) ,
      legend.text=element_text(size=14),
      legend.position="top",
      legend.direction="horizontal")


  if(!is.null(file)){
    pdf(file, width = w, height = h)
    p <- p + geom_text(data=GRSC,aes(label = txt))
    print(p)
    dev.off()
  }
  if(is.null(file)){
    print(p)
  }
}
