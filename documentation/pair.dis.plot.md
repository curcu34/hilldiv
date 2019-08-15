# pair.dis.plot()

The related function pair.dis.plot() uses any of the dissimilarity matrices yielded by pair.dis() (e.g. 1-UqN) to visualize it either as a NMDS chart, a qgraph plot or a heatmap/correlogram.

## Examples
````R
pair.div.q0.L2 <- pair.dis(otu.table[,sort(colnames(otu.table))],qvalue=0,hierarchy=hierarchy.table,level="2")
#Plot NMDS
pair.dis.plot(pair.div.q0.L2$L2_CqN,hierarchy=hierarchy.table,type="NMDS",level=2)
#Plot qgraph
pair.dis.plot(pair.div.q0.L2$L2_CqN,hierarchy=hierarchy.table,type="qgraph",level=2,magnify=TRUE)
````
