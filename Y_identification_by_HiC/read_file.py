def read_file(file):
	if file.endswith('.gz'):
		import gzip
		f = gzip.open(file, 'rt')
		return f
	else:
		f = open(file)
		return f

#ref: https://stackoverflow.com/questions/38655176/reading-in-file-block-by-block-using-specified-delimiter-in-python
def get_groups(f, group_by, sep, idx):
	data = []
	for line in f:
		tmp = line.split(sep)
		if tmp[idx] == group_by:
			if data:
				yield data
				data = []
		data.append(line)
	if data:
		yield data

def store(inFile, keyCol, valCol):
	from collections import defaultdict
	dic = {}
	with open(inFile) as f:
		for line in f:
			line = line.rstrip()
			tmp = line.split('\t')
			dic[tmp[keyCol-1]] = tmp[valCol-1]
	return dic

