for c in 1 3 4 5 7 8 15; do
	python3 exp_module_gh.py c=${c} > ./data/gh/gh_c${c}.csv &
done
