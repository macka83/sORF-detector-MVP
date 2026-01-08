import os
import random
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Konfiguracja Entrez
Entrez.email = "maciek.kala@gmail.com"
# Upewnij się, że Twój klucz API jest poprawny
Entrez.api_key = api_key

def download_genome(accession_id, filename):
    """Pobiera genom z NCBI, jeśli nie istnieje lokalnie."""
    if not os.path.exists(filename):
        print(f"Pobieranie genomu {accession_id}...")
        try:
            with Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text") as handle:
                with open(filename, "w") as f:
                    f.write(handle.read())
            print(f"Zapisano {filename}")
        except Exception as e:
            print(f"Błąd pobierania {accession_id}: {e}")
            return False
    return True

def get_intergenic_orfs(seq_obj, min_aa_len=10, max_aa_len=100):
    """
    Znajduje wszystkie potencjalne ORFy (Start->Stop) na obu niciach DNA.
    Zgodnie z sugestią obsługuje nici forward (+) i reverse (-).
    """
    orfs = []
    
    # Rozszerzona lista start kodonów (dodano CTG)
    start_codons = ['ATG', 'GTG', 'TTG', 'CTG']
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    # Analiza obu nici
    strands_to_check = [
        (seq_obj, "+"),
        (seq_obj.reverse_complement(), "-")
    ]
    
    for current_seq, strand_label in strands_to_check:
        seq_str = str(current_seq)
        for frame in range(3):
            translatable = seq_str[frame:]
            for i in range(0, len(translatable) - 3, 3):
                codon = translatable[i:i+3].upper()
                if codon in start_codons:
                    for j in range(i + 3, len(translatable) - 3, 3):
                        stop_codon = translatable[j:j+3].upper()
                        if stop_codon in stop_codons:
                            orf_dna = translatable[i:j+3]
                            prot_len = len(orf_dna) // 3 - 1
                            if min_aa_len <= prot_len <= max_aa_len:
                                orfs.append((orf_dna, strand_label))
                            break
    return orfs

def extract_negatives_from_genome(gbk_file, strain_name):
    """
    Główna logika ekstrakcji z uwzględnieniem filtrowania promotorów.
    Poprawka: Użycie feature.location.strand zamiast feature.strand [FIX].
    Stosuje asymetryczny masking (2000bp upstream) wg analizy.
    """
    print(f"Przetwarzanie {strain_name}...")
    record = SeqIO.read(gbk_file, "genbank")
    genome_seq = record.seq
    genome_len = len(genome_seq)
    
    occupied = [False] * genome_len
    
    # Parametry maskowania regionów kodujących
    padding_downstream = 30  
    padding_upstream = 2000  
    
    for feature in record.features:
        if feature.type in ["CDS", "rRNA", "tRNA", "gene"]:
            f_start = int(feature.location.start)
            f_end = int(feature.location.end)
            
            # POPRAWKA: strand pobierany z obiektu location
            strand = feature.location.strand 
            
            # Logika upstream zależna od nici
            if strand == -1:
                mask_start = max(0, f_start - padding_downstream)
                mask_end = min(genome_len, f_end + padding_upstream)
            else:
                mask_start = max(0, f_start - padding_upstream)
                mask_end = min(genome_len, f_end + padding_downstream)
                
            for i in range(mask_start, mask_end):
                occupied[i] = True
                
    # Ekstrakcja czystych regionów
    intergenic_seqs = []
    current_chunk = []
    
    for i in range(genome_len):
        if not occupied[i]:
            current_chunk.append(str(genome_seq[i]))
        else:
            if current_chunk:
                intergenic_seqs.append("".join(current_chunk))
                current_chunk = []
    if current_chunk:
        intergenic_seqs.append("".join(current_chunk))
        
    negative_orfs = []
    for region in intergenic_seqs:
        if len(region) >= 30:
            found_orfs = get_intergenic_orfs(Seq(region))
            negative_orfs.extend(found_orfs)
            
    print(f"  -> Znaleziono {len(negative_orfs)} pseudo-sORFów w {strain_name}.")
    return negative_orfs

def validate_shuffled_sequence(seq_str):
    """Odrzuca sekwencje z nienaturalną gęstością STOP kodonów."""
    stop_codons = ['TAA', 'TAG', 'TGA']
    stops = [seq_str[i:i+3] for i in range(0, len(seq_str)-2, 3)]
    stop_count = sum(1 for s in stops if s in stop_codons)
    if stop_count > (len(seq_str) / 3) * 0.15: 
        return False
    return True

def generate_shuffled_decoys(sequences, target_count):
    """Augmentacja danych poprzez tasowanie z walidacją jakości."""
    decoys = []
    current_count = 0
    attempts = 0
    max_attempts = target_count * 10
    
    print(f"Generowanie {target_count} sekwencji syntetycznych z walidacją...")
    
    while current_count < target_count and attempts < max_attempts:
        attempts += 1
        src = random.choice(sequences)
        l = list(src)
        random.shuffle(l)
        shuffled_seq = "".join(l)
        
        if validate_shuffled_sequence(shuffled_seq):
            decoys.append(Seq(shuffled_seq))
            current_count += 1
            
    return decoys

# --- GŁÓWNA KONFIGURACJA ---
TARGET_SIZE = 12000 
genomes = [
    ("U00096.3", "ecoli_k12.gbk"),        # K-12 MG1655, referencyjny lab strain
    ("BA000007.3", "ecoli_o157.gbk"),     # O157:H7 Sakai, patogenny (EHEC)
    ("AE014075.1", "ecoli_cft073.gbk"),   # CFT073, uropatogenny (UPEC)
    ("AP009048.1", "ecoli_w3110.gbk"),    # K-12 W3110, blisko spokrewniony z MG1655
    ("CP001509.1", "ecoli_bl21.gbk"),     # BL21(DE3), do ekspresji białek
    # ("U000606552.1", "ecoli_bw2952.gbk"), # B/W 2952, kliniczny lab strain
    # ("AE005174.1", "ecoli_hs.gbk")        # HS, historyczny reference
]

all_negative_data = [] 

for acc, fname in genomes:
    if download_genome(acc, fname):
        seqs = extract_negatives_from_genome(fname, acc)
        all_negative_data.extend(seqs)

seen = set()
unique_records = []
for seq, strand in all_negative_data:
    if seq not in seen:
        seen.add(seq)
        unique_records.append((seq, strand))

print(f"Łącznie unikalnych biologicznych negatywów: {len(unique_records)}")

final_records = []

for i, (seq_str, strand) in enumerate(unique_records):
    rec = SeqRecord(
        Seq(seq_str),
        id=f"neg_bio_{i}_strand{strand}",
        description="Biological pseudo-sORF"
    )
    final_records.append(rec)

if len(final_records) < TARGET_SIZE:
    missing = TARGET_SIZE - len(final_records)
    biological_seqs_only = [s[0] for s in unique_records]
    shuffled_seqs = generate_shuffled_decoys(biological_seqs_only, missing)
    
    for i, seq_obj in enumerate(shuffled_seqs):
        rec = SeqRecord(
            seq_obj,
            id=f"neg_syn_{i}",
            description="Synthetic shuffled sample (validated)"
        )
        final_records.append(rec)

filename = "large_negative_dataset.fasta"
script_dir = os.path.dirname(os.path.abspath(__file__))
output_file = os.path.join(script_dir, filename)

if final_records:
    count = SeqIO.write(final_records, output_file, "fasta")
    print(f"SUKCES: Zapisano {count} sekwencji do {output_file}.")
