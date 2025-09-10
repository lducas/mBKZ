for c in 1 3 5 15; do
	for seed in {0..4}; do
		python3 exp_slope.py c=${c} ft=dd m=240 > ./data/slopes_240/slope_c${c}_seed${seed}.csv &
		python3 exp_slope.py c=${c} ft=dd m=240 struct=0 > ./data/slopes_240/nostruct_slope_c${c}_seed${seed}.csv &
	done
done
