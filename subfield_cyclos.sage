def is_cyclo(K):
	cc = len(K.roots_of_unity())
	L = CyclotomicField(cc)
	return L.is_isomorphic(K)


def t2(K):
	d = K.degree()
	D = K.disc()
	return (log(abs(D))/(2*d) - log(d)/2)

X = []
for c in range(2, 105):
	if c%4 == 2:
		continue
	K = CyclotomicField(c)
	tK = t2(K)
	# print(t2(K))
	Ls = K.subfields()
	print("conductor %4d = %s subfields:"%(c, factor(c)))
	for L, _, _ in Ls:
		t = t2(L)
		print("  cond %4d, deg %4d : t2 = %.6f"%(c, L.degree(), t), is_cyclo(L))
		if t < 0:
			X.append((t, L.degree(), c, random(), L))
		if t < tK:
			print("subfield has smaller relative discriminant !")
			print("Field: t2=",tK , K)
			print("subfield: t2=", t ,L)
			exit(1)

print("No subfield with smaller relative discriminant encountered.")

print("\n best ones \n")
X.sort()
for (t, d, c, _, L) in X:
	print("original field cond %4d, subfield deg %4d : t2 = %.6f \t"%(c, d, t), "is cyclo:", is_cyclo(L))
