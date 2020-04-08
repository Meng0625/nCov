setwd("C:/Users/Meng/Desktop/nCov/bCellfNet/")

library("igraph")
g <- make_empty_graph(n = 5) %>%
  add_edges(c(1,2, 2,3, 3,4, 4,5)) %>%
  set_edge_attr("color", value = "red") %>%
  add_edges(c(5,1), color = "green")
E(g)[[]]
plot(g)
# ##//['IGHG1', 'IGHA1', 'IGHGP', 'IGHM', 'IGHG3', 'IGHG2','IGHA2', 'IGHD', 'IGHG4', 'None'],
# //[ "#91D1C2", "#DC0000", "#7E6148", "#B09C85"]
# #

color = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")

names(color) <-c('IGHG1', 'IGHA1', 'IGHM', 'IGHG3', 'IGHG2','IGHA2')
par(mfrow=c(3,5),mar=c(0, 0, 3, 0), mgp=c(3, 1, 0), las=0)

# 李博老师的code


for (file in sample_files) {
  data = read.csv(paste0(c("./data",file),collapse = "/"), )
  data = data.frame(data$from,data$to,data$weight)
  graph = graph_from_data_frame(data, directed = FALSE) 
  V(graph)$color <- color[V(graph)]
  E(graph)$width <- 2*data$data.weight
  E(graph)$color <-"grey70"
  E(graph)$frame.color <- "grey30"
  V(graph)$frame.color <- "grey80"
  V(graph)$size <- 85
  V(graph)$label <- NA
  plot(graph, layout=layout_in_circle(graph),margin=0,main=sub('....$','',file))
  # title(paste(), size=10)
  # legend(x=-2.5, y=-.1, c('IGHG1', 'IGHA1', 'IGHM', 'IGHG3', 'IGHG2','IGHA2'), pch=21, col="#777777", pt.bg=color, pt.cex=2, cex=.8, bty="n", ncol=1)
  # sa
}
legend(x=1, y=-.5, c('IGHG1', 'IGHA1', 'IGHM', 'IGHG3', 'IGHG2','IGHA2'), pch=21, col="#777777", pt.bg=color, pt.cex=2, cex=.8, bty="n", ncol=1)