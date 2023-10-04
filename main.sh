inputpath="/home/dengler_yannick/Documents/Isospin_2_analysis/Isospin_2_Analysis/input/isospin_logfiles"

filename="out_scattering_I2"
loglist="./input/isospin_logfiles_list"

find $inputpath  -name $filename  > $loglist
python3 scripts/HDF5.py $loglist
# julia scripts/ensemble_table.jl
# julia scripts/average.jl
# python3 scripts/fitting.py
