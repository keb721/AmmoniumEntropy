#!/bin/bash/

./rungro_3 ${1} 

gmx grompp -f min.mdp -c input_${1}_.gro -p input.top -maxwarn 4 -po mdout_${1}_.mdp -o topol_${1}_.tpr 2>> output_${1}_.log
gmx mdrun -nb cpu -s topol_${1}_.tpr -o traj_${1}_.trr -c co_${1}_.gro -e en_${1}_.edr -g md_${1}_.log 2>> output_${1}_.log
		
thermon=`grep -B1 "Maximum force" output_${1}_.log | head -1 | awk '{print $4}'`
echo ${1} $thermon >> m0.15_gromin_energy_c211_N3_b8.txt 

./mm input_${1}_.gro co_${1}_.gro >> m_0.15_dipangle.txt

rm *_${1}_*
