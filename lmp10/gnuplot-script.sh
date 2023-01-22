echo $1
echo $2
printf "set terminal png\nset output '$1'\nplot '$2','myplot'\n" | gnuplot