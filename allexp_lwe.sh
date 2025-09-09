for seed in {0..14}; do
    python3 exp_lwe.py n=120 verb=0 mBKZ=1 mKannan=1 c=15 ft=dd > ./data/lwe_240/mKan_mBKZ_seed${seed}.csv &
    python3 exp_lwe.py n=120 verb=0 mBKZ=0 mKannan=1 c=15 ft=dd > ./data/lwe_240/mKan_BKZ_seed${seed}.csv &
    python3 exp_lwe.py n=120 verb=0 mBKZ=0 mKannan=0 c=15 ft=dd > ./data/lwe_240/Kan_BKZ_seed${seed}.csv &
done