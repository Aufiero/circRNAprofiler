Version 1.0.0 (Oct 29, 2019)
============================

First release on bioconductor

Version 1.1.2 (development version) (Dic 15, 2019)
==================================================

-   Improved vignettes

-   In plotMotifs function:

    -   To avoid infinity as a value for fold change, 1 was added to
        number of occurences of each motif found in the foreground and
        background target sequences before the normalization.

    -   Added angle param for the rotation angle of the axis labels.
    -   Added removeNegLog2FC param to remove the RBPs having a negative
        log2FC.

-   volcanoPlot function:

    -   Added geneSet param that allows to show only specified host gene
        names.

-   In getMotifs function:

    -   Included possibilities to use motifs of RBPs from MEME database
