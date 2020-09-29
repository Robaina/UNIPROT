## Analyzing the frequency of aminoacids in proteins
Inspired by the work of David McClure (http://dclure.org/) &mdash; where he analyzes the frequency of words across narrative time for thousands of novels, I propose to run a similar study in which the frequencies of aminoacids are analyzed across protein sequence.

We can use for that the entire UniProt (www.uniprot.com) database. Specifically the *Uniparc* section. There, sequences are groupped in taxonomic clusters and separated into fully (human) curated, the *sprot* section, and non-curated, the *trembl* section.

We can download the curated files, they are text files, write a parser in Python or R, and analyze the frequencies across normalized protein length for all organisms and for specific taxons.


## Analogy between protein/DNA sequence evolution and language 
This analogy has been exploited before. However, I was just thinking. What about drawing connections between words in different languages and sequence motifs? Different sequence motifs from different species could have the same evolutionary origin and perform same, similar or different functions. Just like word cognates. But, there are also sequence motifs which do not share an origin and yet perform similar functions. Just like with non-cognate words in different languages. 

What can we ask and try to answer under this framework that has not being asked before? An interestingly universal area of research would be studying the dynamics of meaning adquisition by a sequence space embeded in an interpreting system. Such as a cell or a human collective. In both cases, the sequences alone mean nothing. They are only a bunch of symbols. The only way to tell that these symbols may be part of a system that interpret (and create) them is by looking at deviations from pure randomness. In other words, measuring the (purely sequencial) information that the sequences contain. 

__Note__: Has anyone proposed measurements to quantify information that depends on the system itself? Does this make sense at all? So we would have information solely coming from the sequence and information that appears when the sequence is put in context. Like additional meanings to words? The sequence alone tells us that a word is something that is not random. But a word has different meanings so this is extra information that depennds on the system (i.e. us humans speaking the language).


## Amino acid costs vs frequency 
Cysteine is very costly, however it's very similar to alanine: only one subsitution H -> SH. Incorporating SH is the costly part, since Alanine is rather cheap to produce. Alanine is very frequent while cysteine is not. This is an indication, solely based on cost and frenquence and structure, that the -SH group must be functionally important. And we know it is. Can we apply this analysis to all amino acids?

There must be something to discover/add among all these ideas and concepts. Right now I cannot see it, though. Perhaps another time.