# MicrobiomeViz Tool                                             

<a href="https://drive.google.com/uc?export=view&id=1zETGzCT9MXmVxVBUnbjc2eWJQ5ZDBTPB"><img src="https://drive.google.com/uc?export=view&id=1zETGzCT9MXmVxVBUnbjc2eWJQ5ZDBTPB" style="width: 30%; max-width: 100%; height: auto" title="Microbiome" align="right" />


 Visualize microbiome abundance by uploading tables using R/Shiny.
  
## About the microbiome
  
  A microbiome is a community of micro organisms(bacteria, fungi, viruses) living in a particular ecosystem. 

## Features

- Generate bar plots and pie charts with input from tables
  - *Change treshold level of  with a slider*
 <a href="https://drive.google.com/uc?export=view&id=1u0Ym2pVXIGqQ420Pp7WEXj75Z9GL8lfW"><img src="https://drive.google.com/uc?export=view&id=1u0Ym2pVXIGqQ420Pp7WEXj75Z9GL8lfW" style="width: 100%; max-width: 100%; height: auto" title="rec1" align="left" />
  - *Change type of plot visualization from a dropdown menu (pie chart or bar chart)*
  - *Change taxonomic level from a dropdown menu* 
  - *Panel view with all plots enabled by a checkbox*
- Change between plots from navigation bar menu

## Setup/User manual
1. Go to http://microbiomeviztool.online
2. Upload needed files (*description :bookmark_tabs: is to the right, if new - download 📥 the reference files from the "Read Me" box and upload them instead*)
3. Press the "Add" button
4. Once the plot is visualized, use the controls to the left to filter data and update plots 
5. Add more plots by going to the menu or remove one at a time by pressing the "Remove" button on the plot you want to be deleted

## Running locally (Step by step)
  ##### On Windows:
  1. Open Command Prompt (with administration rights if needed)
  2. Type `ipconfig`

## Future ideas
  - Incorporating other input fomats, such as BIOM

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
