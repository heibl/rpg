#' Plot heatmap of cofidence on MSA
#' @description Plots a heatmap of residue scores values on a MSA
#' @param obj output of guidance, guidance2 or HoT
#' @param file if a path is supplied the alignment is stored in a pdf, if file=NULL (default) alignment is promted
#' @param title title of the plot
#' @param legend.text logical, if TRUE legend is printed
#' @param guidance_score logical, if TRUE a barplot with the GUIDANCE score is printed
#'
#' @author Franz-Sebastian Krah
#' @import ggplot2
#' @import scales
#' @import zoo
#' @importFrom cowplot plot_grid
#' @export

confidence.heatmap <- function(obj, file = NULL, title, legend.text = TRUE,
  guidance_score = TRUE){

  if (!length(grep("score", names(obj)))>0)
    stop("obj is not an output of guidance, guidance2 or HoT")

  txt <- as.vector(as.character(obj$base_msa))
  if (length(grep("scores", names(obj))) > 0)
    obj <- obj$scores

  mat <- data.frame(obj$residue_pair_residue_score, txt)
  rown <- max(mat$residue)
  coln <- max(mat$col)
  res_mat <- matrix(mat$score, nrow = rown, ncol = coln)

  ## scales for PDF
  w <- coln/12
  h <- rown/6

  p <- ggplot(mat,aes(col, residue)) +
    geom_tile(aes(fill = score), colour = "white") +
    scale_fill_gradient2(low = "lightblue", mid = "white",
      high = "purple", midpoint = 0.5, na.value = "gray80") +
    ylab("Sequences (input order)") +
    # xlab("Sites") +
    scale_y_reverse(breaks= pretty_breaks()) +
    ylab("Sequences (input order)")+
    theme_bw()+
    theme(legend.position = "none")


  if(legend.text){
    p <- p + theme(legend.position="top",
      legend.title=element_text(size=10) ,
      legend.text=element_text(size=10),
      legend.direction="horizontal",
      legend.title.align=0.5,
      axis.title.y = element_text(size = rel(1), angle = 90))+
      guides(fill=guide_legend(title.position = "top",
        keywidth = 1, keyheight = 1,
        title=paste("MSA confidence scale", "\n",
          "uncertain <--------> confident", sep="")))
  }

  if(!missing(title)){
    p <- p + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
  }

  if(guidance_score){
    score = obj$residue_pair_column_score
    z <- zoo(score$col_score, score$col)
    z <- merge(zoo(,1:max(score$col)), z)
    score <- data.frame(score = z)
    score <- data.frame(col = 1:dim(score)[1], score)
    score$score[is.na(score$score)] <- 0
    p2 <- ggplot(score, aes(x = col, y = score))+
      geom_bar(stat = "identity", position = position_dodge(),
        fill="darkblue", colour="darkblue")+
      theme_classic()+
      theme(axis.title.y = element_text(size = rel(1), angle = 90))+
      ylab("GUIDANCE score") + xlab("Column")
  }


  if (!is.null(file)){
    p <- p + geom_text(data = mat, aes(label = txt))
    p <- plot_grid(p, p2, ncol = 1, nrow = 2, scale = c(1,1), rel_heights = c(1, 0.3))
    pdf(file, width = w, height = 10)
    print(p)
    dev.off()
  }
  if(is.null(file))
    if(guidance_score)
      p <- plot_grid(p, p2, ncol = 1, nrow = 2, scale = c(1,1),
        rel_heights = c(1, 0.3))
  return(p)
}






