#!/usr/bin/env bash

# author: Donal Burns
# Date: 26/05/2020
###############################################################################
# compile the report
texcount -1 -sum Report.tex > ./Report.sum # word count
sleep 2s
pdflatex Report.tex
pdflatex Report.tex
biber Report
pdflatex Report.tex
pdflatex Report.tex
# remove junk files produced 
sleep 5s # force a pause to ensure the files are removed, without they are left behind for some reason

# remove junk files 
bash ../tex_cleanup.sh ./