if __name__ == '__main__':

    from ParseUniprot import getProteinSequences, getAAFrequencies
    import pickle as pkl

    def saveToPickleFile(python_object, path_to_file='object.pkl'):
        """
        Save python object to pickle file
        """
        out_file = open(path_to_file, 'wb')
        pkl.dump(python_object, out_file)
        out_file.close()

    def readFromPickleFile(path_to_file='object.pkl'):
        """
        Load python object from pickle file.
        Returns python object.
        """
        in_file = open(path_to_file, 'rb')
        python_object = pkl.load(in_file)
        return python_object

    fileNames = ['uniprot_sprot_plants.txt', 'uniprot_sprot_human.txt',
                 'uniprot_sprot_bacteria.txt', 'uniprot_sprot_archaea.txt',
                 'uniprot_sprot_fungi.txt', 'uniprot_sprot_vertebrates.txt',
                 'uniprot_sprot_invertebrates.txt', 'uniprot_sprot_viruses.txt']
    workDir = 'C:\\Users\\tinta\\OneDrive\\Documents\\Projects\\UNIPROT\\'

    amino_acid_frequencies = {}
    for file in fileNames:
        group = file.split('_')[-1].split('.txt')[0]
        print(f'Evaluating {group} database...')
        proteins = getProteinSequences(workDir + file)
        amino_acid_frequencies[group] = getAAFrequencies(proteins, nsections=30)

    saveToPickleFile(amino_acid_frequencies,
                     path_to_file=workDir + 'amino_acid_frequencies.pkl')
