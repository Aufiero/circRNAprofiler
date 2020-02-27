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

Version 1.1.3 (development version) (Dic 19, 2019)
==================================================

-   In plotMotifs function:
    -   Added n param for filtering the motifs.

Version 1.1.4 (development version) (Jen 08, 2020)
==================================================

-   In getMotifs function:
    -   Fixed a bug when trimming the BSJ sequences. E.g. by setting
        width = 6, the BSJ sequences are trimmed so that only 5
        nucleotides are left at each side of the BSJ. If width = 7, then
        6 nucleotides are left at each side of the BSJ etc.

Version 1.1.5 (development version) (Jen 09, 2020)
==================================================

-   In plotMotifs function:
    -   Before only motifs found in the foreground sequences were
        reported. Now also motifs found in the background are reported
        in the output.

Version 1.1.14 (development version) (Jen 09, 2020)
===================================================

-   In mergeBSJunctions function:
    -   Since different detection tools can report sligtly different
        coordinates before grouping the back-spliced junctions, it is
        possible to fix the latter using the gtf file. See param
        fixBSJsWithGTF.
