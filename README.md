# rrapidmarkov

Has this worked?

R code for quickly running and analysing Markov models in health economic evaluation

I wrote this code in mid 2020 after finding the Markov modelling on HEEMOD rather slow.  HEEMOD is a great package, but rather slow for 'vanilla' Markov models.  I finally got round to uploading this here in January 2021.  

Since writing the package, I have became aware of another package on CRAN called 'speedymarkov', which seems to solve exactly the same problems in an almost identical way... and is rather more developed than my version.  As I've used rrapidmarkov in a couple of manuscripts I have uploaded the code here, but from what I've seen speedymarkov seems to be a more complete package so rather than duplicating efforts I'd recommend speedymarkov for general use.  I will continue working on this package as and when though - mainly for my own interest.  Use by anyone else is at your own risk!

I've included userguide and speedcheck vignettes, but R studio is refusing to let me include them in the package for some reason.  I've added the two files as html files in main directory listing above (userguide.html and validation-and-speed-check.html).  To view these, download them and open them with a web browser.  The .rmd files are in the vignettes folder.  (If anyone can shed some light on how to include the vignettes please message me!)  

Ed Wilson

Norwich, UK

Jan 2021

Installation Instructions

To install this package you will need to use the install_github() function in the devtools package:

install.packages("devtools")

library(devtools)

install_github("EdCFWilson/rrapidmarkov")
