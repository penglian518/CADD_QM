#!/bin/bash
#
# This script is used to amplify a single cell to a supercell for the output Trajectory file of CPMD (TRAJEC.xyz)
# Peng Lian 2009-11-1

echo -n "Which Trajectory File to be Amplified? (*.xyz): "
read trj
echo -n "The parameters of PBC (x,y,z): "
read paraPBC
echo $paraPBC | grep "," > /dev/null 2>&1
if [ "$?" -ne "0" ]; then
 echo "Parameters of PBC must be spaced by \",\"!"
 echo -n "The parameters of PBC (x,y,z): "
 read paraPBC
 echo $paraPBC | grep "," > /dev/null 2>&1
 if [ "$?" -ne "0" ]; then
  echo "Parameters of PBC must be spaced by \",\"!"
  exit
 fi
fi
a=`echo $paraPBC | gawk -F, '{print $1}'`
b=`echo $paraPBC | gawk -F, '{print $2}'`
c=`echo $paraPBC | gawk -F, '{print $3}'`

echo -n "How many cells will be placed along X axis? "
read na
echo -n "How many cells will be placed along Y axis? "
read nb
echo -n "How many cells will be placed along Z axis? "
read nc

#trj=last10000.xyz
#a=5.1832
#b=6.2357
#c=8.5181
#na=3
#nb=3
#nc=3


Tmp=/tmp/rdf_cal$$
mkdir $Tmp

cat -n $trj | grep STEP > $Tmp/noriginsteps.tmp
gawk '{print $1}' $Tmp/noriginsteps.tmp > $Tmp/forcount.tmp
gawk '{print $2"	    "$3}' $Tmp/noriginsteps.tmp > $Tmp/originsteps.tmp
nstep=`cat $Tmp/forcount.tmp | wc -l`
n1=`head -1 $Tmp/forcount.tmp`
n2=`head -2 $Tmp/forcount.tmp | tail -1`
deltan=`expr $n2 - $n1`
natom=`expr $deltan - 2`
newnatom=`echo "$na*$nb*$nc*$natom" | bc`
head -$deltan $trj | tail -$natom | gawk '{print $1}' > $Tmp/labels.tmp

#echo $nstep $n1 $n2 $deltan $natom

s=1
while [ "$s" -le "$nstep" ]; do
 nhead=`echo "$s*$deltan" | bc`
 ntail=$natom
 head -$nhead $trj | tail -$ntail > $Tmp/astep.tmp
 
 fna=1
 while [ "$fna" -le "$na" ]; do
  deltax=`echo "scale=6; $fna*$a" | bc`
  gawk '{printf "%10.6f\n", $2+"'$deltax'" }' $Tmp/astep.tmp > $Tmp/astep.x
   fnb=1
   while [ "$fnb" -le "$nb" ]; do
    deltay=`echo "scale=6; $fnb*$b" | bc`
    gawk '{printf "%10.6f\n", $3+"'$deltay'" }' $Tmp/astep.tmp > $Tmp/astep.y
     fnc=1
     while [ "$fnc" -le "$nc" ]; do
      deltaz=`echo "scale=6; $fnc*$c" | bc`
      gawk '{printf "%10.6f\n", $4+"'$deltaz'" }' $Tmp/astep.tmp > $Tmp/astep.z
      paste $Tmp/labels.tmp $Tmp/astep.x $Tmp/astep.y $Tmp/astep.z >> $Tmp/newtrj.tmp
      
      fnc=`expr $fnc + 1`
     done
    fnb=`expr $fnb + 1`
  done
  fna=`expr $fna + 1`
 done

echo "          $newnatom" >> $Tmp/newtrj.xyz
head -$s $Tmp/originsteps.tmp | tail -1 >> $Tmp/newtrj.xyz
cat $Tmp/newtrj.tmp >> $Tmp/newtrj.xyz

rm -fr $Tmp/newtrj.tmp
echo "Amplifying for step $s"
s=`expr $s + 1`
done
mv $Tmp/newtrj.xyz ./NEW_TRAJEC_amplified.xyz
rm -fr $Tmp
