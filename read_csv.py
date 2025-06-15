
def read_csv(filename):
	f = open(filename, "r")
	lines = f.read().splitlines()
	headers = lines[0].replace("\t", "").replace(" ", "").split(",")
	c = len(headers)
	data = [[] for h in range(c)]	
	for line in lines[1:]:
		l = line.split(",")
		if len(l) != c:
			continue
		for i in range(c):
			try:
				data[i].append(eval(l[i]))
			except:
				data[i].append(None)
	return {headers[i]:data[i] for i in range(c)}

