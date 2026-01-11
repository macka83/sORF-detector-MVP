# %%
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import gffutils  # pip install gffutils
from BCBio import GFF
from dotenv import dotenv_values

config = dotenv_values(".env")
# Konfiguracja Entrez
Entrez.email = "maciek.kala@gmail.com"
# Upewnij się, że Twój klucz API jest poprawny
Entrez.api_key = config["api_key"]

# Pobierz handle FASTA
handle = Entrez.efetch(db="nuccore", id="U00096.3", rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

# Zapisz do pliku
SeqIO.write(record, "ecoli_k12.fasta", "fasta")

print(f"Zapisano: {len(record.seq)} nt")  # ~4641652 bp
handle = Entrez.efetch(db="nuccore", id="U00096.3", rettype="gb", retmode="text")
with open("ecoli.gbk", "w") as f:
    f.write(handle.read())


in_file = "ecoli.gbk"  # GenBank z U00096.3
out_file = "ecoli.gff3"

with open(in_file) as in_handle, open(out_file, "w") as out_handle:
    GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)

# 1. Pobierz i załaduj genom + anotacje
genome = list(SeqIO.parse("ecoli_k12.fasta", "fasta"))[0]
db = gffutils.create_db(
    "ecoli.gff3",
    ":memory:",
    keep_order=True,
    merge_strategy="create_unique",
    force=True,
)

# 2. Intergeniczne regiony
intergenic = []
genes = db.features_of_type("gene")
for inter in db.interfeatures(
    genes, new_featuretype="intergenic"
):  # lub custom subtract genes
    seq = genome.seq[inter.start - 1 : inter.stop]
    if len(seq) > 100:  # filtr
        intergenic.append(str(seq))

# 3. Generuj negatywy (shuffling z GC match, bez Shine-Dalgarno)
negatives = []
gc = 0.5  # E. coli avg
probs = np.array([(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])  # A/T bias
for i, parent in enumerate(intergenic[:10000]):  # duża ilość
    for _ in range(10):  # wiele na region
        length = np.random.randint(30, 101) * 3  # nt, wielokrotność 3
        shuffled = "".join(np.random.choice(["A", "C", "G", "T"], size=length, p=probs))
        if "AGGAGG" not in shuffled[-20:]:  # brak RBS upstream
            negatives.append(SeqRecord(Seq(shuffled), id=f"neg_{i}_{_}"))

SeqIO.write(negatives, "ecoli_negatives.fasta", "fasta")
