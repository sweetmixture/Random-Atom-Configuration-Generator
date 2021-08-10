#!/bin/python
mx=7000							# Number of randomised cluster user wants
rm -rf top_structures RACG_statistics log		# Delete expected output files, if they exist in the directory
mkdir top_structures					# Output .xyz files in ~
touch RACG_statistics					# structure_code - energy summary: this file data may be sorted on demend (in order to remove duplicated structures)
							# can be don by using 'remove_dup.py'. 
echo "total_trial $mx" > log				

for ((i=1; i<=$mx; ));	do

	python MAIN_RACG.py > rand_res.out		# record temporal racg output (at the step)
	
	if_success=$( grep "RESULT_TAG:" rand_res.out | awk '{print $2}')			# check if random config generation & pre-optimisation succeeded
	
	if [ $if_success == "SUCCESS" ]; then							# if success keyword foudn
		energy=$(grep "SCF DONE" *.xyz | awk '{print $3}')				# save the config
		name=B"$i".xyz
		mv *.xyz $name
		printf "%12.12s\t%12.6f\n" $name $energy >> RACG_statistics			# record structure_code - energy
		mv *.xyz top_structures
		((i++))				# GIVE INCREMENT WHEN ONLY SUCCESS
		sort -k 2 RACG_statistics > tmp
		mv tmp RACG_statistics
	fi
done
