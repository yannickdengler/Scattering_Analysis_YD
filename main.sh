inputpath="/home/dengler_yannick/Documents/Isospin_2_analysis/Isospin_2_Analysis/input/isospin_logfiles"

filename="out_scattering_I2"
loglist="./input/isospin_logfiles_list"

# find $inputpath  -name $filename  > $loglist
python3 dev/HDF5.py $loglist
python3 dev/basic_analysis.py
python3 dev/energy_levels.py
# python3 dev/basic_analysis.py
