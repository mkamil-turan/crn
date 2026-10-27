from crn.fetch import Gencode

gc = Gencode()
#gc.update()

#Bir geni almak için gerekli sysmbol kodunu yaz
g1 = gc.gene("MEFV")
#g1 = gc.gene("TBP")
#g1 = gc.gene("RNA18S1")

#Tüm transkripler ve onlara ait id nuımaraları
print(g1.transcripts.keys())

#Bunlardan birini al
#tid = "ENST00000956137.1
#tid = 'ENST00001089022.1'
tid = 'ENST00000339854.9'

#buna ait rna kodunu bu şekilde al
rna = g1.transcripts[tid]

# CDS 'yi çekebilrisn
s, f = rna.features.get('cds')[0]
mRNA = rna[s-1:f-3]
print(mRNA)

proteinler = [p for p in g1.proteins.values()
              if p.features.first("transcript_id") == tid]
print(proteinler[0])

print(mRNA.translation())

# Ya da doğrudan: transkript ile proteini eşleyen CDS
print(rna.match_orf(proteinler[0]))
