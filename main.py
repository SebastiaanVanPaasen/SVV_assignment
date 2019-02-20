def read_table(filename):  # tables must be in form: 1 row header (description), separation by tabulators
    file = open(filename, "r+")
    file = file.readlines()
    n_colums = len(file[0].split('\t'))
    data=np.array(ndmin=2)
    del file[0]
    for line in file:
        l = line.split('\t')


    for i in (ps1, ps2, t):
        del i[0]
    for i in range(len(t)):
        floatt.append(float(t[i]))
        floatps1.append(float(ps1[i]))
        floatps2.append(float(ps2[i]))
    return (floatt, floatps1, floatps2)
