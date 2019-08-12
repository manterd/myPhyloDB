set_panel_size <- function(p=NULL, g=ggplotGrob(p), file=NULL,
                           margin = unit(1,"mm"),
                           width=unit(4, "cm"),
                           height=unit(4, "cm")){

  panels <- grep("panel", g$layout$name)
  panel_index_w<- unique(g$layout$l[panels])
  panel_index_h<- unique(g$layout$t[panels])
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)

if(getRversion() < "3.3.0"){

   # the following conversion is necessary
   # because there is no `[<-`.unit method
   # so promoting to unit.list allows standard list indexing
   g$widths <- grid:::unit.list(g$widths)
   g$heights <- grid:::unit.list(g$heights)

   g$widths[panel_index_w] <-  rep(list(width),  nw)
   g$heights[panel_index_h] <- rep(list(height), nh)

} else {

   g$widths[panel_index_w] <-  rep(width,  nw)
   g$heights[panel_index_h] <- rep(height, nh)

}

  if(!is.null(file))
    ggsave(file, g,
           width = convertWidth(sum(g$widths) + margin,
                                unitTo = "in", valueOnly = TRUE),
           height = convertHeight(sum(g$heights) + margin,
                                  unitTo = "in", valueOnly = TRUE))

  invisible(g)
}

print.fixed <- function(x) grid.draw(x)


make_biom <- function(data, sample_metadata=NULL, observation_metadata=NULL, id=NULL, matrix_element_type="int"){
  # The observations / features / OTUs / rows "meta" data table
  if(!is.null(observation_metadata)){
    rows = mapply(list, SIMPLIFY=FALSE, id=as.list(rownames(data)),
                  metadata=alply(as.matrix(observation_metadata), 1, .expand=FALSE, .dims=TRUE))
  } else {
    rows = mapply(list, id=as.list(rownames(data)), metadata=NA, SIMPLIFY=FALSE)
  }
  # The samples / sites / columns "meta" data table
  if(!is.null(sample_metadata)){
    columns = mapply(list, SIMPLIFY=FALSE, id=as.list(colnames(data)),
                     metadata=alply(as.matrix(sample_metadata), 1, .expand=FALSE, .dims=TRUE))
  } else {
    columns = mapply(list, id=as.list(colnames(data)), metadata=NA, SIMPLIFY=FALSE)
  }
  # Convert the contingency table to a list
  datalist = as.list(as.data.frame(as(t(data), "matrix")))
  names(datalist) <- NULL
  # Define the list, instantiate as biom-format, and return
  # (Might eventually expose some of these list elements as function arguments)
  format_url = "http://biom-format.org/documentation/format_versions/biom-1.0.html"
  return(biom(list(id=id,
                   format = "Biological Observation Matrix 1.0.0-dev",
                   format_url = format_url,
                   type = "OTU table",
                   generated_by = sprintf("biom %s", packageVersion("biom")),
                   date = as.character(Sys.time()),
                   matrix_type = "dense",
                   matrix_element_type = matrix_element_type,
                   shape = dim(data),
                   rows = rows,
                   columns = columns,
                   data = datalist)
  ))
}