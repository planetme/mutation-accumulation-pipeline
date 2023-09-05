import glob as glob

pairsDict = dict()

for sampleNumber in range(28, 78+1):
    reads = glob.glob(("/mnt/researchdrive/nsharp2/c_albicans_MA/sequence_data/experimental_dataset/" + "*S" + str(sampleNumber) + "*" + ".fastq.gz"))
    pairsDict[sampleNumber] = reads

print(pairsDict)