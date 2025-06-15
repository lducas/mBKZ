for c in 1 3 4 5 7 8 15 16; do
	python3 exp_skewness_index.py c=${c} > ./data/skewness_index/skewness_index_c${c}.csv &
done
