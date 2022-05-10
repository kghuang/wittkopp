import csv
csv_file = open('network_out.tsv')
csv_reader = csv.reader(csv_file, delimiter='\t')
line = 0 # corresponds to column
counts = [] # keeps track of counts for each gene
titles = [] # names of the genes

for row in csv_reader:
    # for first row, record all gene names
    if line == 0: 
        line += 1
        titles = row[1:-19]
        # initialize counts for each gene
        for j in range(len(titles)):
            counts.append({"c": 0, "n": 0})
        continue
    i = 0
    for col in row:
        if i == 0: # skip first column
            i = 1
            continue
        if i >= 2446: # skip the last 19 columns
            break
        # count calculations (i - 1 because we skip first row above)
        counts[i + 1]["c"] += (row[-1] == "Yes" and col == "1")
        counts[i + 1]["n"] += int(col)
        i += 1
    line += 1


results = [titles, []]
i = 0
# calculate proportions
for t in results[0]:
    if counts[i]["n"] == 0:
        results[1].append(0)
    else:
        results[1].append(float(counts[i]["c"]) / float(counts[i]["n"]))
    i += 1
# make new file with results
out_file = open("results.csv", "w")
for row in results:
    for item in row:
        out_file.write(str(item) + ",")
    out_file.write("\n")
