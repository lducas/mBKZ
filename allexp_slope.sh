for seed in {0..20}; do
	python3 exp_lwe.py n=120 ft=dd verb=0 struct=1 > ./data/lwe_120/struct_seed${seed}.csv &
	python3 exp_lwe.py n=120 ft=dd verb=0 struct=0 > ./data/lwe_120/nostruct_seed${seed}.csv &
done
