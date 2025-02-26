This is a template directory for projects using R targets and Seurat v5 to analyze 10X data.



## Metadata

Metadata for each sample should be provided in data/metadata.csv. This should be a comma-separated .csv-file (i.e. US format, not European). The table can contain any number of columns, but should at least have one column, "Sample", where the values match the names of the data files to be read (e.g. sample1 for sample1.h5 or the directory sample1/).

### Metadata for doublet identification

If doublets should be predicted using DoubletFinder and removed, the metadata table should also contain a column labeled "expectedCells", indicating the number of cells expected based on the number of cells loaded into each GEM well. If multiple samples were loaded into the same GEM well and later demultiplexed using hashtags, the column "gemWell" should also be present, with unique values for each GEM well. Note that the values are not important, as long as they are unique. As an exmple, if samples A and B were loaded into the same GEM well, in proportions expected to give 6 000 and 4 000 cells, respectively, and sample C was prepared in another well with an expected output of 10 000 cells, the metadata table could be as follows:

| Sample | expectedCells | gemWell |
|---|---|---|
| SampleA | 6000 | FirstWell |
| SampleB | 4000 | FirstWell |
| SampleC | 10000 | SecondWell |
