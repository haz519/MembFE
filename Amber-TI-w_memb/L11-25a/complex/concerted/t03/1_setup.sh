#!/bin/sh

#windows=$(seq 0.000, 0.091, 0.182, 0.273)

: > mbar.dat

for w in 0.0000 0.0479 0.1151 0.2063 0.3161 0.4374 0.5626 0.6839 0.7937 0.8850 0.9521 1.0000; do
  if [ \! -x $w ]; then
    mkdir $w
  fi

  echo $w >> mbar.dat

  sed -e "s/%L1%/$w/" min.in.tmpl > $w/min.in.tmpl

  sed -e "s/%L1%/$w/" eq.in.tmpl > $w/eq.in.tmpl
  sed -e "s/%L1%/$w/" ti.in.tmpl > $w/ti.in.tmpl
  (cd $w; ln -sf ../mbar.dat)
done
