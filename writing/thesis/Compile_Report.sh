#!/usr/bin/env bash

# author: Donal Burns
# Date: 26/05/2020
###############################################################################
# compile the report
texcount -1 -sum Donal_Burns_CMEE_MSc_2020.tex > ./Donal_Burns_CMEE_MSc_2020.sum # word count
sleep 2s
pdflatex Donal_Burns_CMEE_MSc_2020.tex
pdflatex Donal_Burns_CMEE_MSc_2020.tex
biber Donal_Burns_CMEE_MSc_2020
pdflatex Donal_Burns_CMEE_MSc_2020.tex
pdflatex Donal_Burns_CMEE_MSc_2020.tex
# remove junk files produced 
sleep 5s # force a pause to ensure the files are removed, without they are left behind for some reason

# remove junk files 
bash ../tex_cleanup.sh ./