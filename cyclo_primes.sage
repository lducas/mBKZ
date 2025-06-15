firsts = 10000
print("cyclo_primes = dict([])")
for c in [1, 3, 4, 5, 7, 8, 9, 15, 16, 27]:
	L = []
	K = CyclotomicField(c)
	p = 1
	while len(L) < firsts or L[firsts-1] > p:
		p = next_prime(p)
		Ps = K.primes_above(p)
		for P in Ps:
			L.append(P.norm())
		L.sort()
		L = L[:firsts]
	print("cyclo_primes[%d] = %s"%(c, str(L)))
