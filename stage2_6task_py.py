# -*- coding: utf-8 -*-
"""stage2_6task.py

Original file is located at
    https://colab.research.google.com/drive/1cw_BBdPJljVWps8nkMrH2Y_YPN0qAtbu
"""

# Transcriptomics task on RNA sequencing

# import RNA dataset
import pandas as pd

rna_dataset = ("https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt")
rna_df = pd.read_csv(rna_dataset, sep="\s+")
print(rna_df)
# padj is the adjusted p-value, comparing gene expression between two groups there's a chance of getting false positives. these p-values are adjusted to control the false positives, so padj represents the p-value after these adjustment
# foldchange is the ratio of expression levels between the two groups
# log2foldchange represents the magnitude of change in gene expression between two groups. it is calculated as the logarithm base2 of the fold change



# generate a volcano plot
# a volcano plot is a type of scatter plot used to quickly identify genes, proteins, or other variables that exhibit large magnitude changes that are also statistically significant. They help to quickly identify the most interesting features in a large dataset
# for a volcano plot x-axis:log2foldchange yaxis:-log10(padj)
import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px
import numpy as np

# calculate -log10(padj)
rna_df['-log10(padj)'] = -np.log10(rna_df['padj'])
print(rna_df)

# generate a volcano plot
fig = px.scatter(rna_df, x='log2FoldChange', y='-log10(padj)', hover_name='Gene', title = "Volcano Plot of Gene Expression")
fig

# assuming the "rna_df" contains the entire dataset
# determine the upregulated gene
# the rna_df was defined previoously and needs to be accessible in this cell
upregulated_genes = rna_df[(rna_df['log2FoldChange'] > 1) & (rna_df['pvalue'] < 0.01)]
print(upregulated_genes)

# determine the downregulated gene
downregulated_genes = rna_df[(rna_df['log2FoldChange'] < -1) & (rna_df['pvalue'] < 0.01)]
print(downregulated_genes)

# functions of the top five upregulated gene according to GeneCards
top_upregulated_genes = upregulated_genes.nlargest(5, 'log2FoldChange')
print(top_upregulated_genes)
# DTHD1 gene: encodes a protein which contains a death domain, it functions in signalling pathway and formation of signalling complexes
# EMILIN2 gene: involved in angiogenesis, and positive regulation of platelet aggregation
# PI16 gene: predicted to be involve in negative regulation of peptidase activity
# C4orf45 gene: encodes the acidic form of complement factor 4, part of the classical activation pathway
# FAM180B gene: predicted to be located in extracellular region

# functions of the top five downregulated gene according to GeneCards
top_downregulated_genes = downregulated_genes.nlargest(5, 'log2FoldChange')
print(top_downregulated_genes)

# FAM46B: negative regulation of apoptotic process
# HS3ST3A1: key components in generating a myriad of distinct heparan sulfate fine structures that carry out multiple biologic activities
# PMEL: encodes a protein enriched in melanosomes, plays an essential role in the structural organisation of premelanosomes
# TNFAIP6: encodes a secretory protein that contains a hyaluronan-binding domain which is involved in extracellular matrix stability and cell migration
# COL4A2: encodes one of the six subunits of type IV collagen, the major structural components of basement membranes
