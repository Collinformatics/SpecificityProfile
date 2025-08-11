import random


# Original sequences
inputSeqExp = "GTTATTTTACAACTTGAACGTGTT"
inputSeqBg =  "AAGACTATCACCAATTGGTATCAG"

# Number of variants
N=500


def generateVariants(sequence, mutationOdds=4, numVariants=50):
    seq5Prime = 'AAAGGCAGT'
    seq3Prime = 'GGTGGAAGT'
    bases = ['A', 'T', 'C', 'G']
    variants = [('variant_0', f'{seq5Prime}{sequence}{seq3Prime}')]

    while len(variants) < numVariants:
        var = list(sequence)
        for index in range(0, len(sequence), 3):
            codon = sequence[index:index + 3]
            if random.randint(1, 10) >= mutationOdds:
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
variantsExp = generateVariants(inputSeqExp, 4, numVariants=N)
variantsBg = generateVariants(inputSeqBg, 9, numVariants=N)

# Save to FASTA and FASTQ
pathExp = "data/variantsExp.fastq"
saveSeqs(variantsExp,  pathExp)
pathBg = "data/variantsBg.fasta"
saveSeqs(variantsBg,  pathBg)

print(f"Saved {N} variants at:\n"
      f"     {pathExp}\n"
      f"     {pathBg}")
