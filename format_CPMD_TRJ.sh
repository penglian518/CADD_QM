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

A=`echo "scale=6; 1000000 * $a" | bc | gawk -F. '{print $1}'`
B=`echo "scale=6; 1000000 * $b" | bc | gawk -F. '{print $1}'`
C=`echo "scale=6; 1000000 * $c" | bc | gawk -F. '{print $1}'`

amin=`echo "scale=0; -1 * $A/2" | bc`
amax=`echo "scale=0; 1 * $A/2" | bc`
bmin=`echo "scale=0; -1 * $B/2" | bc`
bmax=`echo "scale=0; 1 * $B/2" | bc`
cmin=`echo "scale=0; -1 * $C/2" | bc`
cmax=`echo "scale=0; 1 * $C/2" | bc`


#trj=last10000.xyz
#a=5.1832
#b=6.2357
#c=8.5181
#na=3
#nb=3
#nc=3


Tmp=/tmp/cpmd_trj_format$$
mkdir $Tmp

cat -n $trj | grep STEP > $Tmp/noriginsteps.tmp
gawk '{print $1}' $Tmp/noriginsteps.tmp > $Tmp/forcount.tmp
gawk '{print $2"	    "$3}' $Tmp/noriginsteps.tmp > $Tmp/originsteps.tmp
nstep=`cat $Tmp/forcount.tmp | wc -l`

n1=`head -1 $Tmp/forcount.tmp`
n2=`head -2 $Tmp/forcount.tmp | tail -1`
deltan=`expr $n2 - $n1`		#how many lines each step will take in the trj file

natom=`expr $deltan - 2`	#number of atoms
#newnatom=`echo "$na*$nb*$nc*$natom" | bc`

head -$deltan $trj | tail -$natom | gawk '{print $1}' > $Tmp/labels.tmp

#echo $nstep $n1 $n2 $deltan $natom

s=1
while [ "$s" -le "$nstep" ]; do
 nhead=`echo "$s*$deltan" | bc`
 ntail=$natom
 head -$nhead $trj | tail -$ntail > $Tmp/astep.tmp
 
 gawk '{printf "%10.6f\n", 1000000*$2}' $Tmp/astep.tmp | gawk -F. '{print $1}' > $Tmp/astep.x
 gawk '{printf "%10.6f\n", 1000000*$3}' $Tmp/astep.tmp | gawk -F. '{print $1}' > $Tmp/astep.y
 gawk '{printf "%10.6f\n", 1000000*$4}' $Tmp/astep.tmp | gawk -F. '{print $1}' > $Tmp/astep.z

 for num in `cat $Tmp/astep.x`
 do
  if `echo $num | grep - > /dev/null 2>&1`; then
   while [ "$num" -lt "$amin" ]; do
    num=`expr $num + $A`
   done
  else
   while [ "$num" -gt "$amax" ]; do
    num=`expr $num - $A`
   done
  fi
  echo "scale=6; $num/1000000" | bc >> $Tmp/newstep.x
 done

 for num in `cat $Tmp/astep.y`
 do
  if `echo $num | grep - > /dev/null 2>&1`; then
   while [ "$num" -lt "$bmin" ]; do
    num=`expr $num + $B`
   done
  else
   while [ "$num" -gt "$bmax" ]; do
    num=`expr $num - $B`
   done
  fi
  echo "scale=6; $num/1000000" | bc >> $Tmp/newstep.y
 done

 for num in `cat $Tmp/astep.z`
 do
  if `echo $num | grep - > /dev/null 2>&1`; then
   while [ "$num" -lt "$cmin" ]; do
    num=`expr $num + $C`
   done
  else
   while [ "$num" -gt "$cmax" ]; do
    num=`expr $num - $C`
   done
  fi
  echo "scale=6; $num/1000000" | bc >> $Tmp/newstep.z
 done
 
 gawk '{printf "%10s \n", $1}' $Tmp/newstep.x > $Tmp/newstep.x.gawk
 gawk '{printf "%10s \n", $1}' $Tmp/newstep.y > $Tmp/newstep.y.gawk
 gawk '{printf "%10s \n", $1}' $Tmp/newstep.z > $Tmp/newstep.z.gawk

 paste $Tmp/labels.tmp $Tmp/newstep.x.gawk $Tmp/newstep.y.gawk $Tmp/newstep.z.gawk >> $Tmp/newtrj.tmp

 echo "          $natom" >> $Tmp/newtrj.xyz
 head -$s $Tmp/originsteps.tmp | tail -1 >> $Tmp/newtrj.xyz
 cat $Tmp/newtrj.tmp >> $Tmp/newtrj.xyz

 rm -fr $Tmp/newtrj.tmp
 rm -fr $Tmp/newstep.x*
 rm -fr $Tmp/newstep.y*
 rm -fr $Tmp/newstep.z*
 echo "formatting for step $s"
 s=`expr $s + 1`
done
mv $Tmp/newtrj.xyz ./NEW_TRAJEC_format.xyz
rm -fr $Tmp
