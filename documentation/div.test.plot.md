# div.test.plot()

Visual comparison between the diversity levels of two or multiple groups of samples, with the option to show the significance values of the pairwise posthoc analyses.

## Examples
````R
contrast.div.q0 <- div.test(otu.table,qvalue=0,hierarchy=hierarchy.table,posthoc=TRUE)
colours <- c("#35a849","#9d1923","#f7ab1b","#ed7125","#cc4323","#b6d134","#fcee21","#085ba7")

#Plots without posthoc analysis results
div.test.plot(contrast.div.q0)
div.test.plot(contrast.div.q0,chart="box")
div.test.plot(contrast.div.q0,chart="jitter",colour=colours)
div.test.plot(contrast.div.q0,chart="violin",colour=c("#35a849","#9d1923","#f7ab1b","#ed7125","#cc4323","#b6d134","#fcee21","#085ba7"))

#Plots with posthoc analysis results
div.test.plot(contrast.div.q0,chart="box",colour=colours,posthoc=TRUE)
div.test.plot(contrast.div.q0,chart="jitter",colour=colours,posthoc=TRUE,threshold=0.5)
````
<img align=left src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.test.plot.png" width="400" title="div.test.plot() with pairwise comparisons">
<img src="https://github.com/anttonalberdi/DiverHill/blob/master/figures/div.test.plot.violin.png" width="400" title="Violin div.test.plot() with pairwise comparisons">
