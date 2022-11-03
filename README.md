# MicrobiomeViz Tool                                             

<a href="https://drive.google.com/uc?export=view&id=1zETGzCT9MXmVxVBUnbjc2eWJQ5ZDBTPB"><img src="https://drive.google.com/uc?export=view&id=1zETGzCT9MXmVxVBUnbjc2eWJQ5ZDBTPB" style="width: 30%; max-width: 100%; height: auto" title="Click to enlarge picture" align="right" />


 Visualize plots based on a human microbiome by uploading tables using R/Shiny.
  
## About the microbiome
  
  A microbiome is a community of micro organisms living in an ecosystem. This project will be focusing on the human one (the micro organisms living in us).

## Features

- Generate multiple plots based on information from tables
- Filter data, reset plots, visualize more by using sliders, dropdown and checkbox controls: 
  - *Change treshold level with a slider*
  - *Change type of plot visualization from a dropdown menu (pie chart or bar chart)*
  - *Change taxonomic level from a dropdown menu* 
  - *Panel view with every plot by checking and unchecking a checkbox*
- Change between plots from navigation bar menu
- Remove any unwanted tab with a plot 

## Setup/User manual
1. Go to the website
2. Upload needed files (*description :bookmark_tabs: is to the right, if new - download 📥 the reference files from the "Read Me" box and upload them instead*)
3. Press the "Add" button
4. Once the plot is visualized, use the controls to the left to filter data and reload plots 
5. Add more plots by going to the menu or remove one at a time by pressing the "Remove" button on the plot you want to be deleted

## Running locally (Step by step)
  
## Future ideas
  - App functioning with BIOM files

## How are the tables supposed to look like?
  ###### Metadata file: 
  |  | sample_id | disease_stat | ... |
  | --- | --- | --- | --- |
  | 1. | ID of the sample | Type of the case | ... |
  
  ###### Subsample.shared(OTU) file: 
  | Group | Otu0001 | ... |
  | --- | --- | --- |
  | DA00006 | value | ... |
  
  ###### Taxonomy file: 
  | OTU | Taxonomy | ... |
  | --- | --- | --- |
  | Otu0001 | Bacteria(100);Bacteroidetes(100);... | ... |
