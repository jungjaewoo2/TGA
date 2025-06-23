#!/bin/sh 
./ascii.sh $1 > jnk.sgml
sgml2txt -man jnk.sgml
mv jnk.man $1
