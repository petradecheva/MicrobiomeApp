# MicrobiomeViz Tool                                             

<a href="https://drive.google.com/uc?export=view&id=1zETGzCT9MXmVxVBUnbjc2eWJQ5ZDBTPB"><img src="https://drive.google.com/uc?export=view&id=1zETGzCT9MXmVxVBUnbjc2eWJQ5ZDBTPB" style="width: 30%; max-width: 100%; height: auto" title="Microbiome" align="right" />


 Visualize microbiome abundance by uploading tables using R/Shiny.
  
## About the microbiome
  
  A microbiome is a community of micro organisms (bacteria, fungi, viruses) living in a particular ecosystem. 
## Features

- Generate bar plots and pie charts with input from tables:
  - *Change treshold level (of what?) with a slider (for the purpose of what?)*
 <a href="https://drive.google.com/uc?export=view&id=1u0Ym2pVXIGqQ420Pp7WEXj75Z9GL8lfW"><img src="https://drive.google.com/uc?export=view&id=1u0Ym2pVXIGqQ420Pp7WEXj75Z9GL8lfW" style="width: 100%; max-width: 100%; height: auto" title="rec1" align="left" />
  - *Change type of plot visualization from a dropdown menu (pie chart or bar chart)*
  - *Change taxonomic level from a dropdown menu* 
  - *Panel view with all plots enabled by a checkbox*
- Change between plots from navigation bar menu

## Setup/User manual
1. Go to the website
2. Upload needed files (explanation is to the right in the app)
3. Press "Add" 
4. Once the plot is visualized, use the controls to the left to filter and update 
5. Add more plots by going to the menu or remove one at a time by pressing the "Remove" button on the plot you want to be deleted

## Running locally (Step by step)
  
## Future ideas
  - Incorporating other input formats, such as BIOM 

## Input format:

  ###### Metadata file: 
  |  | sample_id | disease_stat | ... |
  | --- | --- | --- | --- |
  | 1. | Sample ID | Sample type | ... |
  
  ###### Subsample.shared(OTU) file: 
  | label | Group | numOTUs | Otu0001... | ... |
  | --- | --- | --- | --- | --- |
  | Label | Sample type | OTU's sum | number of specific OTU | ... |
  
  ###### Taxonomy file: 
  | OTU | Size | Taxonomy | ... |
  | --- | --- | --- | --- |
  | OTU label | Number of sequences in the OTU| Consensus taxonomy | ... |
