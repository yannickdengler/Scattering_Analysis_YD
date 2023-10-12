inputpath="/home/dengler_yannick/Documents/Isospin_2_analysis/Isospin_2_Analysis/input/isospin_logfiles"

filename="out_scattering_I2"
loglist="./input/isospin_logfiles_list"

find $inputpath  -name $filename  > $loglist
python3 dev/HDF5.py $loglist
echo "HDF5 done!"
python3 dev/basic_analysis.py
echo "basic_analysis done!"
python3 dev/energy_levels.py
echo "energy_levels done!"
python3 dev/infinite_volume.py
echo "infinite_volume done!"
python3 dev/phase_shift.py
echo "phase_shift done!"