import os

from crn.fetch import NcbiSession, Ncbi

# API anahtarı koda yazılmaz: NCBI_API_KEY ortam değişkeninden okunur (yoksa anonim, 3 istek/sn).
session = NcbiSession(api_key=os.environ.get("NCBI_API_KEY"))

ncbi = Ncbi(session=session, organism="human")
mefv = ncbi.gene("4210")          # ya da ncbi.gene("MEFV")
print(mefv)

for tid, rna in mefv.transcripts.items():
    start, end = rna.features.first("cds")       # ör. (47, 2392)
    cds = rna[start - 1:end]                     # kodlayan bölge
    print(tid, len(rna), rna.features.first("protein_id"))


prt = mefv.protein("NP_000234.1")                # 781 aa
print(prt.features.first("transcript_id"))       # NM_000243.3
print(mefv.features.first("map_location"))