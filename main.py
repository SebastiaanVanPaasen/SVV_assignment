def read_table(filename):  # tables must be in form: 1 row header (description), separation by tabulators
    file = open(filename, "r+")
    file = file.readlines()
    n_colums = len(file[0].split('\t'))
    del file[0]
    n_rows = len(file)
    data=np.zeros(shape=(n_colums,n_rows))
    for i in range(n_rows):
        l = file[i].split('\t')
        for j in range(n_colums):
            data[j,i]=l[j]
    return data
