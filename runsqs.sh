#!/bin/bash

rm -rf Ti-*

# Define variables
for i in 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.20;do
	j=$(printf "%.0f" $(echo "100 * $i" | bc -l))
	mkdir Ti-${j}pctV
	cd Ti-${j}pctV
	variable1=3.23741703254859
	variable2=${i}
	variable3=$(echo "1 - $variable2" | bc)  # Calculating 1 - variable2
	replication=5 # Should be divisible by solute percentage, 5x5x5 replication of 2 atoms unitcell makes 250 atoms, which is divisible by 0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.20
	first_NN=$(echo "(sqrt(3)/2) * $variable1" | bc -l)
	second_NN=$(echo "1 * $variable1" | bc -l)
	random_distance=$(awk -v min="$first_NN" -v max="$second_NN" 'BEGIN{srand(); print min+rand()*(max-min)}')
    tolerance=0.01

	# Create rndstr.in file : primitive cell
	cat <<EOF > rndstr.in
$variable1 $variable1 $variable1 90 90 90
-.5 .5 .5
.5 -.5 .5
.5 .5 -.5
.000000 .000000 .000000 Ti=$variable3,V=$variable2
EOF

	# Create sqscell.out file : SuperCell

	cat <<EOF > sqscell.out 
1

$replication 0 0
0 $replication 0
0 0 $replication
EOF

	corrdump -l=rndstr.in -ro -noe -nop -clus -2=${random_distance}
	mcsqs -rc -tol=${tolerance}

	# Complie C++ code to convert sqs to POSCAR
	cp ../sqs2poscar.cpp .
	cp ../poscar2data.py .
	c++ ./sqs2poscar.cpp -o ./sqs2poscar

	for input_file in bestsqs-*.out; do
  	if [ -e "$input_file" ]; then
    	./sqs2poscar "$input_file"
  	else
    	echo "File $input_file does not exist."
  	fi
	done
	python3 poscar2data.py
	cd ..
done
