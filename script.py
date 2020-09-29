if __name__ == '__main__':
    
    from ParseUniprot import plotfreqAAs
    # plt.style.use('dark_background')
    FileName = 'uniprot_sprot_plants.txt'
    WorkDir = 'C:\\Users\\tinta\\OneDrive\\Documents\\Projects\\UNIPROT\\'
    sol = plotfreqAAs(WorkDir+FileName, nsections=30)
