violintheme = function(){
  theme(axis.text = element_text(size = 8),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 8))
}


addColors = function(plot, scale_type = "colour", layer = NULL){
  if(!exists("get_color_list")){
    source("https://raw.githubusercontent.com/jenfransson/r_resources/refs/heads/main/code/color_list.R")
  }
    if(is.null(layer)){
      if(! scale_type %in% names(plot$mapping))
        stop(paste0(scale_type," is not mapped in the given plot"))
      map_col = 
        plot$mapping[[scale_type]][[2]]
      map_data = plot$data
    }else{
      if(! scale_type %in% names(plot$layers[[layer]]$mapping))
        stop(paste0(scale_type," is not mapped in the given plot"))
      map_col = 
        plot$layers[[layer]]$mapping[[scale_type]][[2]]
      if(class(plot$layers[[layer]]$data) != "waiver"){
        map_data = plot$layers[[layer]]$data
      }else{
        map_data = plot$data
      }
    }
  
  if(class(map_col) == "call"){
    map_col = map_col[[3]]
  }
  
  nvalues = length(unique(map_data[[map_col]]))
  plot & get(paste0("scale_",scale_type,"_manual"))(
    values = get_color_list(nvalues))
}
