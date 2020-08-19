#!/usr/bin/env bash

# author: Donal Burns
# Date: 26/05/2020
###############################################################################
# compile the report
texcount -1 -sum Burns_Donal_CMEE_MSc_2020.tex > ./Burns_Donal_CMEE_MSc_2020.sum # word count
sleep 2s
pdflatex Burns_Donal_CMEE_MSc_2020.tex
pdflatex Burns_Donal_CMEE_MSc_2020.tex
biber Burns_Donal_CMEE_MSc_2020.tex
pdflatex Burns_Donal_CMEE_MSc_2020.tex
pdflatex Burns_Donal_CMEE_MSc_2020.tex
# remove junk files produced 
sleep 5s # force a pause to ensure the files are removed, without they are left behind for some reason

# remove junk files 
bash ../tex_cleanup.sh ./