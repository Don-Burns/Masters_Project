#!/bin/bash
# to remove the "trash" files generated when compiling latex `.tex` files

rm $1*.aux
rm $1*.bbl
rm $1*.blg
rm $1*.log
rm $1*.out
rm $1*.toc
rm $1*.dvi
rm $1*.gz
rm $1*.fls
rm $1*.fdb_latexmk
rm $1*.synctex.gz
rm $1*.bcf
rm $1*.run.xml