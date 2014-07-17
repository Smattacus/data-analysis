function A = read_PMTdata(filename)
f1 = fopen(filename);
A = fread(f1, '*uint');