#!/bin/bash
#

export LD_LIBRARY_PATH=$AMBERHOME/lib_gnu:$LD_LIBRARY_PATH

mbarstates=$(cat mbar.dat | wc -l)
mbarlambda=$(tr '\n' ' ' < mbar.dat)
jidx=3
curr=ti$(printf "%03i" $jidx)
prev=$(($jidx - 1))
prev=ti$(printf "%03i" $prev)

if [ $jidx -eq 1 ]; then
  sed -e "s/%L2%/$mbarstates/" -e "s/%L3%/$mbarlambda/" min.in.tmpl > $curr.min.in
  sed -e "s/%L2%/$mbarstates/" -e "s/%L3%/$mbarlambda/" eq.in.tmpl > $curr.eq.in
  sed -e "s/%L2%/$mbarstates/" -e "s/%L3%/$mbarlambda/" ti.in.tmpl > $curr.prod.in
  # minimization
  $AMBERHOME/bin/pmemd.cuda \
    -i $curr.min.in -p ../concerted.parm7 -c ../concerted.rst7 \
    -O -o $curr.min.out -inf $curr.min.info -e $curr.min.en -r $curr.min.rst7 \
    -x $curr.min.nc -l $curr.min.log -AllowSmallBox

  # equilibration, heating
  $AMBERHOME/bin/pmemd.cuda \
    -i $curr.eq.in -p ../concerted.parm7 -c $curr.min.rst7 \
    -O -o $curr.eq.out -inf $curr.eq.info -e $curr.eq.en -r $curr.eq.rst7 \
    -x $curr.eq.nc -l $curr.eq.log -AllowSmallBox

  # first TI run
  $AMBERHOME/bin/pmemd.cuda \
    -i $curr.prod.in -p ../concerted.parm7 -c $curr.eq.rst7 \
    -O -o $curr.prod.out -inf $curr.prod.info -e $curr.prod.en -r $curr.prod.rst7 \
    -x $curr.prod.nc -l $curr.prod.log -AllowSmallBox
  wait

else
  sed -e "s/%L2%/$mbarstates/" -e "s/%L3%/$mbarlambda/" ti.in.tmpl > $curr.prod.in
  # next TI run
  $AMBERHOME/bin/pmemd.cuda \
    -i $curr.prod.in -p ../concerted.parm7 -c $prev.prod.rst7 \
    -O -o $curr.prod.out -inf $curr.prod.info -e $curr.prod.en -r $curr.prod.rst7 \
    -x $curr.prod.nc -l $curr.prod.log -AllowSmallBox
  wait

fi

