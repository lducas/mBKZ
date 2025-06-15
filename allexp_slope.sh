for c in 1 3 4 5 8 15 16; do
	for seed in {0..4}; do
		python3 exp_slope.py c=${c} ft=dd m=240 > ./data/slopes_240/slope_c${c}_seed${seed}.csv &
	done
done
