# AI_NET_BIO
This pipeline is dedicated to biomarker discovery, starting from a raw gene expression matrix, passing through normalization, harmonization, DEGs cross-validation, and building an adjacency network and a random forest model for hub disease prediction.
Please take a look at the complete illustrated pipeline below
<img width="1536" height="2784" alt="AI_NET_BIO" src="https://github.com/user-attachments/assets/fb65963e-5a02-4050-91fc-1e8d4ccad0da" />


## 🚀 What codes are included in the project
* **Complete Pipeline for Harmonization with plot-packed exploratory analysis**: harmonization.R file includes R codes that can be used to harmonize gene expression matrices from microarray experiments from any platform (Affymatrix, Illumina, Agilent)up to DEGs identification using limma.
* **2500-iteration Crossvalidation**: Crossvalidation+Overlap_ Calculation.R is code to crossvlaidate DEGs by running 2500 iterations and calculating the overlap coefficient, then visualizing the overlap between iterations.
* **Correlation and Network Construction**: cal_correlation+adjacency_network_buiding.R is a simple code to address the correlation of your cross-validated DEGs expression and construct a co-expression adjacency matrix ready to be analyzed by the MONET package.
* **Biological relevance through functional predication**: enrichment_modules_disease_genes.R is a pathway and GO enrichment code exhibited from cluster profiles, that includes various visualization plots for enriched pathway including circular plot generation.
* **The AI integration**: build_RF_model.R is an R-based code to build a random forest model to get the top-ranked hubs from your filtered modules of the network produced at the previous stage. Also, it can be run independently, with code ready to plot the ROC curve to assess the model's performance.

## 📦 Installation

This is a complete R-based pipeline, with commented instructions for clear illustration of usage.

## ## 📄 Paper
This project is a part of the paper entitled: Identifying behavior regulatory leverage over mental disorders transcriptomic network hubs toward lifestyle-dependent psychiatric drugs repurposing
(https://link.springer.com/article/10.1186/s40246-025-00733-w)
We'd love you to get use of the pipeline and enjoy coding ^_^.
Please cite it if you are going to use the code.
