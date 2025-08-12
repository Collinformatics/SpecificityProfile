import random


# Original sequences
seqDNA = 'GTTATTTTACAACTTGAACGTGTT' # Starting protein sequence
seq5Prime = 'AAAGGCAGT' # 5' flanking sequence
seq3Prime = 'GGTGGAAGT' # 3' flanking sequence

# Dataset parameters

N=1 # Number of variants
oddMutationExp = 0 # Percent chance of mutating a codon in experimental set
oddMutationBg = 0

setInit = False
if setInit:
    pathExp = "data/variantsExp.fastq"
    pathBg = "data/variantsBg.fasta"
else:
    seqDNA = 'GCTTTAATTCAAATTGATAATGCT'
    pathExp = "data/variantsExp2.fastq"
    pathBg = "data/variantsBg2.fasta"


def generateVariants(sequence, mutationOdds=4, numVariants=50):

    bases = ['A', 'T', 'C', 'G']
    variants = [('variant_0', f'{seq5Prime}{sequence}{seq3Prime}')]

    while len(variants) < numVariants:
        var = list(sequence)
        for index in range(0, len(sequence), 3):
            codon = sequence[index:index + 3]
            if  random.randint(0, 100) <= mutationOdds:
                bp = codon[random.randint(0, 2)]
                bpNew = random.choice([b for b in bases if b != bp])
                var[index] = bpNew
        if var != sequence:
            name = f"variant_{len(variants)}"
            subCassette = seq5Prime + "".join(var) + seq3Prime
            variants.append((name, subCassette))

    return variants


def saveSeqs(variants, fileName):
    if '.fastq' in fileName:
        with open(fileName, 'w') as fastq:
            for name, seq in variants:
                quality = ''.join(chr(random.randint(53, 73)) for _ in seq)
                fastq.write(f"@{name}\n{seq}\n+\n{quality}\n")
    else:
        with open(fileName, 'w') as fasta:
            for name, seq in variants:
                fasta.write(f">{name}\n{seq}\n")


# Generate variants
variantsExp = generateVariants(seqDNA, mutationOdds=oddMutationExp, numVariants=N)
variantsBg = generateVariants(seqDNA, mutationOdds=oddMutationBg, numVariants=N)

# Save to FASTA and FASTQ
saveSeqs(variantsExp,  pathExp)
saveSeqs(variantsBg,  pathBg)

print(f"Saved {N} variants at:\n"
      f"     {pathExp}\n"
      f"     {pathBg}")
