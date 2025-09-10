for seed in {0..29}; do
    python3 exp_lwe.py n=136 verb=0 mBKZ=1 mKannan=1 c=15 ft=dd > ./data/lwe_272/mKan_mBKZ_seed${seed}.csv &
    python3 exp_lwe.py n=136 verb=0 mBKZ=0 mKannan=1 c=15 ft=dd > ./data/lwe_272/mKan_BKZ_seed${seed}.csv &
    python3 exp_lwe.py n=136 verb=0 mBKZ=0 mKannan=0 c=15 ft=dd > ./data/lwe_272/Kan_BKZ_seed${seed}.csv &
done