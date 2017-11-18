#!/bin/bash
rm result
touch result
if [ $# -eq 0 ]
	then
		echo "No arguments supplied"
else
	if [ "$1" -gt "0" ]
		then 
			echo $1

		EXT=cfg

		for i in configurations/*; do
		    if [ "${i}" != "${i%.${EXT}}" ];then
		        echo "Executing with : $i" >> result
		        for j in `seq 1 1 "$1"`;
					do
						echo "$j"
				    	sequentiel/genetic data/graphe4.txt ${i} >> result
					done

		        

		    fi
		done
	fi	
fi


