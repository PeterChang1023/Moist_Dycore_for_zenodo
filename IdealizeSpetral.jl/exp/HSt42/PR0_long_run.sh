#!bin/sh
"""
This script is designed for long-duration simulations, producing output results at intervals specified by 'space_day'. It helps prevent excessive memory usage by periodically saving results instead of waiting until the end of the simulation.
"""

start_day=0
final_day=20000
space_day=25

for i in `seq $start_day $space_day $final_day`
do
L=0
echo $i"day"
echo -n $L > Latent_heat.txt


if [ $i -eq 0 ]; then
	rm -rf HSt42_${L}
	mkdir HSt42_${L}

    rm -rf warmstart_cp_HSt42_${L}
    mkdir warmstart_cp_HSt42_${L}
	echo -n $space_day > HSt42_${L}/day_interval.txt
	echo -n "None" > HSt42_${L}/firstday_file.txt
else
	echo -n $space_day > HSt42_${L}/day_interval.txt
	echo -n "warmstart_${L}.dat" > HSt42_${L}/firstday_file.txt
fi

julia Run_HS.jl
L=0 

if [ -f "HSt42_${L}/warmstart_${L}.dat" ] && [ $i -lt $final_day ]; then
	cp "HSt42_"${L}"/warmstart_${L}.dat" "warmstart_${L}.dat"
	cp "HSt42_"${L}"/warmstart_${L}.dat" "warmstart_cp_HSt42_${L}/warmstart_${L}_$(($i+$space_day))th_day.dat"
	mv "HSt42_"${L}"/all_L"${L}".dat" "HSt42_"${L}"/RH80_PR"$L"_"$final_day"day_startfrom_"$i"day_final.dat"
	echo 'warmstart file exists.'
elif [ -f "HSt42_${L}/warmstart_${L}.dat" ] && [ $i -eq $final_day ]; then
	mv "HSt42_"${L}"/all_L"${L}".dat" "HSt42_"${L}"/RH80_PR"$L"_"$final_day"day_startfrom_"$i"day_final.dat"
	cp "HSt42_"${L}"/warmstart_${L}.dat" "warmstart_cp_HSt42_${L}/warmstart_${L}_$(($i+$space_day))th_day.dat"
	mv "HSt42_"${L}"/warmstart_${L}.dat" "HSt42_"${L}"/HSt42_"${L}"RH80_PR"$L"_"$final_day"day_startfrom_"$i"day_final.dat"
	echo "All files have completed!!!"
fi
done
