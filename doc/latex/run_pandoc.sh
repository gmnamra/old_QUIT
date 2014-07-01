#!/bin/bash

pandoc -s --toc --mathml --highlight-style monochrome DESPOT.tex -o ../html/despot.html --to=html5 --bibliography=/Users/Tobias/Documents/Kings/Bibliography/rat_bib.bib --css=/Users/Tobias/Code/QUIT/man/DESPOT/test.css
