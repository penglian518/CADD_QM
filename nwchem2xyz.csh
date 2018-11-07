#!/bin/csh -f

if ($#argv < 1) then
  echo Usage:
  echo 'nwchem2xyz.csh nwchem_output_file.out [all]'
  exit
endif

set dollar = \$
set xyz = $1:r.xyz

if (-e $xyz) rm $xyz

#set lcoords = `grep -n 'Output coordinates' $1 | \
#    sed -e 's/Output coordinates in angstroms (scale by  1.889725989 to convert to a.u.)//g' -e 's/:/ /g'`
set lcoords = `grep -n '          Step ' $1 | sed -e 's/Step....//g' -e 's/:/ /'`

if ($2 == 'all') then
  set i = 0
else
  @ i = $#lcoords - 1
endif

while ( $i < $#lcoords )
  @ i++
  @ j = $i - 1
#  @ lprint = $lcoords[$i] + 4
  @ lprint = $lcoords[$i] + 11 
  sed -n "$lprint, /^ *$dollar/ p" $1 | sed '/^ *$/d' > tmp 
  set atoms = `wc -l tmp`
  echo $atoms[1] >> $xyz
  echo Step $j >> $xyz
  cat tmp | awk '{printf "%-5s %14.9f %14.9f %14.9f\n", $2, $4, $5, $6}' >> $xyz 
  rm tmp
end
if ($2 == 'all') then
  echo wrote $i snapshots
else
  echo wrote last snapshot only. Use \'all\' argument for trajectory.
endif


