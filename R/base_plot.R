

create_plot_data <- function(counts,name,coord)
{
  dt = data.table(x=coord,y =counts,condition = name)
  return(dt)
}
