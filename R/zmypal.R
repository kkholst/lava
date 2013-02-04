mypal <- function(set=TRUE,...) {
  col <- c("black","darkblue","darkred","goldenrod","mediumpurple",
  "seagreen","aquamarine3","violetred1","salmon1",
  "lightgoldenrod1","darkorange2","firebrick1","violetred1", "gold")
  if (!set) return(col)
  palette(col)
}
                         
