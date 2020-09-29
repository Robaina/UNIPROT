#FileName = 'uniprot_sprot_archaea.txt'
#WorkDir = 'C:\\Users\\tinta\\OneDrive\\Documents\\Projects\\Frenquency_of_aminoacids\\'
# TODO: Would be interesting to compute all protein costs per each protein in E. coli and
# plot the distribution. Also, perhaps we can predict the distribution of costs using 
# techniques from statistical physics. We would have to constraint the average to some
# value, the physiological limit of protein synthesis. 
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import re
import warnings

def getProteinSequences(DataFile):

    # Read protein sequence data file and extract protein sequences

    # Read master protein data file
    with open(DataFile) as file:
        data = file.read()

    # Isolate protein sequences (using key: CRC64 with 12 empty spaces before start)
    pstart = [m.start() for m in re.finditer('CRC64', data)]
    pend = [m.start() for m in re.finditer("\n//", data)]

    proteins = {}
    for idx in range(len(pstart)):
        proteins[idx] = re.sub(' ','',re.sub('\n','',data[pstart[idx]+12:pend[idx]]))
    
    return proteins


def plotfreqAAs(DataFile, nsections=4):

    # Read protein sequence data file and compute the frequency of each aminoacid
    # accross the protein sequence and all proteins. It also plots the resulting
    # bar plots for the 20 aminoacids.

    proteins = getProteinSequences(DataFile)

    aminoacids = {'G':'Glycine','A':'Alanine','L':'Leucine','M':'Methionine',
                  'F':'Phenylalanine','W':'Tryptophan','K':'Lysine','Q':'Glutamine',
                  'E':'Glutamic acid','S':'Serine','P':'Proline','V':'Valine',
                  'I':'Isoleucine','C':'Cysteine','Y':'Tyrosine','H':'Histidine',
                  'R':'Arginine','N':'Asparagine','D':'Aspartic acid','T':'Threonine'}
    aas = list(aminoacids.keys())
    naas = len(aas)

    # Compute the frequency of each aminoacid across protein sequence
    aminoacid_freqs = {}
    for aa in aas:
        aa_counts = np.zeros(nsections)
        for k in range(len(proteins)):
            if len(proteins[k]) > 0:
                subprots = np.array_split(list(proteins[k]),nsections)
                for idx in range(len(subprots)):
                    aa_counts[idx] += subprots[idx].tolist().count(aa)
        aminoacid_freqs[aa] = 100*aa_counts/aa_counts.sum()

    # Plot histograms for each aminoacids
    fig = plt.figure(figsize = (20, 10))
    plt.suptitle(DataFile.split('_')[-1].split('.txt')[0].capitalize(),fontsize = 20)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
    fig.text(0.5, 0.04, 'normalized distance to sequence start', va='center', ha='center', fontsize=14)
    fig.text(0.02, 0.5, 'frequency', va='center', ha='center', rotation='vertical', fontsize=14)
    for idx,aa in enumerate(aas):
        plt.subplot(4,5,idx+1)
        plt.bar(range(nsections), aminoacid_freqs[aa])
#         plt.ylabel('frequency (%)', fontsize = 10)
        plt.title(aminoacids[aa], fontsize = 12)
    plt.show()

    return({'aa_freqs':aminoacid_freqs,'protein_seqs':proteins})

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def plot_aminoacid_freq(protein):

    # Computes the histogram depicting the frequency of every aminoacid in a
    # protein. The argument 'protein' is a list containing the one-letter code
    # of its amminoacids.

    
    aminoacids = {'G':'Glycine','A':'Alanine','L':'Leucine','M':'Methionine',
                  'F':'Phenylalanine','W':'Tryptophan','K':'Lysine','Q':'Glutamine',
                  'E':'Glutamic acid','S':'Serine','P':'Proline','V':'Valine',
                  'I':'Isoleucine','C':'Cysteine','Y':'Tyrosine','H':'Histidine',
                  'R':'Arginine','N':'Asparagine','D':'Aspartic acid','T':'Threonine'}

    aa_counts = []
    aas = list(aminoacids.keys())
    protein_seq = list(protein)
    protein_len = len(protein_seq)
    for aa in aas:
        aa_counts.append(100*protein_seq.count(aa)/protein_len)

    sorted_idx = np.argsort(-np.array(aa_counts))
    aa_counts.sort(reverse = True)
    aa_names = [aminoacids[aas[i]] for i in sorted_idx]
    naas = 20

    # Plot the results in bar plot
    plt.figure(figsize = (12,6))
    plt.bar(range(naas), aa_counts)
    plt.xticks(range(naas),aa_names, rotation = 75, fontsize = 18)
    plt.ylabel('frequency (%)', fontsize = 18)
    plt.title('Protein length: ' + str(len(protein_seq)) + ' aas', fontsize = 20)
    plt.show()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def getProtSeq(DataFile, genes):

    # This function retrieves the protein sequence corresponding to the genes in the
    # the argument. Genes is a list containing Escherichia coli gene indentifiers of
    # the form 'bxxx'.

    # Read master protein data file
    with open(DataFile) as file:
        data = file.read()

    # Find genes
    ProtSeqs = {}
    for gene in genes:
        featureIdx = [m.start() for m in re.finditer('KEGG; eco:' + gene, data)]
        IdxStart = data[featureIdx[0]:featureIdx[0]+10000].find('CRC64')
        IdxEnd = data[featureIdx[0]:featureIdx[0]+10000].find('\n//')
        ProtSeqs[gene] = list(re.sub(' ','',re.sub('\n','',data[featureIdx[0]+IdxStart+12:featureIdx[0]+IdxEnd])))

    return(ProtSeqs)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def getProteinCost(proteinSequence):

    # Average number of high energy phosphate groups required to synthetize one
    # molecule, according to Akashi & Gojobori, 2002, PNAS. ProtSeq is a list
    # containing the aminoacids (single letter code) in the protein.

    #aa_costs = {'G':11.7,'A':11.7,'S':11.7,'D':12.7,'N':14.7,'E':15.3,'Q':16.3,
    #            'T':18.7,'P':20.3,'V':23.3,'C':24.7,'L':27.3,'R':27.3,'K':30.3,
    #            'I':32.3,'M':34.3,'H':38.3,'Y':50,'F':52,'W':74.3}

    # Costs as computed in Krick et al. 2014 (the above corrected by decay rate
    # units of Phosphate/Time)

    aa_costs = {'G':12,'A':12,'V':47,'L':55,'P':61,'I':65,'S':70,'E':77,
            'R':109,'T':112,'D':114,'Q':130,'N':147,'F':208,'K':242,
            'Y':350,'M':446,'H':536,'C':741,'W':892}
                
    return(sum([aa_costs[aa] for aa in proteinSequence if aa in aa_costs.keys()]))


def getAllProteinsCosts(DataFile):
    """
    Get protein costs for all proteins in the data file, together with their 
    sequence length
    """
    warnings.filterwarnings('ignore')
    proteinCosts, proteinSizes = [], []
    proteins = getProteinSequences(DataFile)
    for protein in proteins.values():
        psize = len(protein)
        proteinSizes.append(psize)
        proteinCosts.append(psize * getProteinCost(protein))
     
    cost = np.log(np.array(proteinCosts))
    proteinSize = np.array(proteinSizes)
    
    plt.figure(figsize=(12, 9))
    plt.xlabel('log(cost)')
    plt.ylabel('counts')
#     plt.hist(np.log(np.array(proteinCosts)))
    sns.distplot(cost, hist=True, kde=True,
                         color='blue', kde_kws={'shade': True, 'linewidth': 3},
                         hist_kws={'edgecolor': 'black'})
    plt.show()
    
    plt.figure(figsize=(12, 9))
    plt.xlabel('size')
    plt.ylabel('counts')
    sns.distplot(proteinSize, hist=True, kde=True,
                         color='green', kde_kws={'shade': True, 'linewidth': 3},
                         hist_kws={'edgecolor': 'black'})
    plt.show()
    
    warnings.resetwarnings()
    return proteinCosts
    
        
   
