#!/bin/sh 
./ascii.sh $1 > jnk.sgml
sgml2latex -d jnk.sgml
mv jnk.dvi $1.dvi
