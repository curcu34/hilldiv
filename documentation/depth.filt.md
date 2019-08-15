# depth.filt()

Filter samples based on a minimum sequencing depth.

## Arguments
| Arguments | Description |
| ------------- | ------------- |
| otutable | An OTU table (matrix/data.frame) indicating the absolute OTU abundances of multiple samples. Columns must refer to samples and rows to OTUs.
| threshold | A number indicating the minimum sequencing depth required to keep the sample. |

## Examples
````R
depth.filt(otu.table,5000)
depth.filt(otu.table,threshold=20000)
````
