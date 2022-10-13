# tippi



## License

This work is licensed under a Creative Commons Attribution 4.0 International License

http://creativecommons.org/licenses/by/4.0/

![http://creativecommons.org/licenses/by/4.0/](https://i.creativecommons.org/l/by/4.0/88x31.png)



## Structure & analysis order

The repo is structured as:

```
.
├── dataprep
│   ├── graphs
│   └── outdata
├── inference
│   ├── data
│   ├── graphs
│   └── outdata
└── model
    ├── data
    ├── graphs
    ├── indata
    └── outdata
```

Analyses should be run in this order:

- *dataprep* : preparation of pre/post data for inference. Also contains some utility functions.
- *inference* : statistical analyses of pre/post data - can all be run and post-processed using the bash script.
- *model* : the main modelling analysis
  - *modeldata.R* : prepares effect estimates from the inference analysis and makes parts of Table 1
  - *modeloutcomes.R* : main results for each sensitivity scenario
  - *runmodel.sh* : bash script to run all modelling analyses at once (on multicore machine with plenty memory)

*resultsupload.R* in the top level parses and uploads results to googlesheets for inclusion in the article. The uploads require authetication (authors only).

### Required packages ###

The following R packages are required:

- inference: rstanarm
- plotting: ggplot2, scales, ggpubr, ggthemes
- modelling: HEdtree, discly (install using `devtools::install_github('petedodd/packagename')` )
- data manipulation: data.table
- utilities: glue, here,
- data I/O: googlesheets4, readxl
