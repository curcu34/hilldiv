# div.test.plot()

Visual comparison between the diversity levels of two or multiple groups of samples. The 'chart' argument enables selecting between boxplot, violin plot and jitter plot. The 'stat' argument enables plotting pairwise mean comparisons across groups.

## Examples
Visual comparison between the diversity levels of two or multiple groups of samples. The 'chart' argument enables selecting between boxplot, violin plot and jitter plot. The 'stat' argument enables plotting pairwise mean comparisons across groups.
````R
contrast.div.q0 <- div.test(otu.table,qvalue=0,hierarchy=hierarchy.table)
colours <- c("#35a849","#9d1923","#f7ab1b","#ed7125","#cc4323","#b6d134","#fcee21","#085ba7")

#Box plot
div.test.plot(contrast.div.q0)
div.test.plot(contrast.div.q0,chart="box")

#Jitter plot
div.test.plot(contrast.div.q0,chart="jitter",colour=colours)

#Violin plot
div.test.plot(contrast.div.q0,chart="violin",colour=c("#35a849","#9d1923","#f7ab1b","#ed7125","#cc4323","#b6d134","#fcee21","#085ba7"))

#Pairwise mean comparison statistical significances
div.test.plot(contrast.div.q0,chart="jitter",stat=TRUE,flip=TRUE)
div.test.plot(contrast.div.q0,stat=TRUE,comb=list(c("Myotis myotis","Myotis capaccinii")),symbol=TRUE)
````
<img align=left src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.test.plot.png" width="400" title="div.test.plot() with pairwise comparisons">
<img src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.test.plot.violin.png" width="400" title="Violin div.test.plot() with pairwise comparisons">
