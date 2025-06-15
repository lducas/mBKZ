for c in 1 3 4 5 8 15 16; do
	for seed in {0..4}; do
		python3 exp_profile.py c=${c} ft=dd beta=64 m=160 > ./data/profiles_d160_beta64/prof_c${c}_seed${seed}.csv &
	done
done
