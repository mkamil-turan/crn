from crn.io import DNA, Strand


dna1 = DNA(seq="atgcatgcatgc", strand=Strand.FivePrime)
dna1.name = "Teo1"
dna1.features.add("gene","teo1")
dna1.features.add("type","coding")
print(dna1)
print(type(dna1))
dna1.to_crn("test/dna1.crn")


print("Buda crn dosyasından okudun")
dna2 = DNA.from_crn("test/dna1.crn")
print(dna2)
print(type(dna2))
