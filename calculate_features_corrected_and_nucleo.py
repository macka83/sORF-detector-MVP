from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import numpy as np
import os
import warnings
from collections import Counter

# Ignorujemy ostrzeżenia o nieznanych aminokwasach (np. X)
warnings.filterwarnings("ignore")

def calculate_aliphatic_index(seq):
    """
    Oblicza indeks alifatyczny białka.
    Wzór: X(Ala) + 2.9 * X(Val) + 3.9 * (X(Ile) + X(Leu))
    """
    L = len(seq)
    if L == 0: return 0
    a_count = seq.count('A')
    v_count = seq.count('V')
    il_count = seq.count('I') + seq.count('L')
    
    return ((a_count + 2.9 * v_count + 3.9 * il_count) / L) * 100

def calculate_gc_content(seq):
    """Oblicza zawartość GC w sekwencji DNA."""
    if not seq: return 0
    seq_u = str(seq).upper()
    gc = seq_u.count('G') + seq_u.count('C')
    return gc / len(seq)

def calculate_hexamer_entropy(seq):
    """Oblicza entropię Shannona dla heksamerów (struktura sekwencji)."""
    seq = str(seq).upper()
    if len(seq) < 6: return 0
    
    hexamers = [seq[i:i+6] for i in range(len(seq)-5)]
    if not hexamers: return 0
    
    counter = Counter(hexamers)
    total = len(hexamers)
    freqs = np.array(list(counter.values())) / total
    
    entropy = -np.sum(freqs * np.log2(freqs + 1e-10))
    # Normalizacja przez max możliwą entropię dla znalezionych unikalnych heksamerów
    return entropy / np.log2(len(counter) + 1) if len(counter) > 0 else 0

def calculate_codon_bias(seq):
    """
    Uproszczony bias kodonów: AT-richness na 3. pozycji (wobble).
    Wiele bakteryjnych sORF wykazuje specyficzny bias w tej pozycji.
    """
    seq = str(seq).upper()
    if len(seq) < 3: return 0
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    if not codons: return 0
    
    wobble_at = sum(1 for c in codons if c[-1] in ['A', 'T'])
    return wobble_at / len(codons)

def calculate_physicochemical_features(fasta_file):
    results = []
    print(f"Przetwarzanie sekwencji z pliku: {fasta_file}...")
    
    rejected_len = 0
    rejected_err = 0
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        # 1. Translacja DNA -> Białko (table 11 = Bacterial)
        try:
            protein_seq = record.seq.translate(table=11, to_stop=True)
        except Exception:
            continue
        
        # 2. Filtracja długości (sORF: 10-100 aminokwasów)
        if len(protein_seq) < 10 or len(protein_seq) > 100:
            rejected_len += 1
            continue
            
        str_protein = str(protein_seq).replace("X", "").replace("J", "").replace("Z", "").replace("*", "")
        str_dna = str(record.seq)
        
        if len(str_protein) < 10:
            continue

        try:
            analysed_seq = ProteinAnalysis(str_protein)
            
            # 3. Budowa słownika cech
            features = {
                "id": record.id,
                # Cechy DNA (NOWE - Zgodnie z rekomendacją)
                "dna_sequence": str_dna,
                "dna_len": len(str_dna),
                "orf_gc": calculate_gc_content(str_dna),
                "orf_entropy": calculate_hexamer_entropy(str_dna),
                "orf_codon_bias": calculate_codon_bias(str_dna),
                
                # Placeholder na kontekst (wymagane przez sorf_level1_update.py)
                "upstream_seq": "", 
                "downstream_seq": "",
                
                # Cechy Białkowe (Istniejące - Poprawne)
                "prot_len": len(str_protein),
                "molecular_weight": analysed_seq.molecular_weight(),
                "gravy": analysed_seq.gravy(),
                "aromaticity": analysed_seq.aromaticity(),
                "instability": analysed_seq.instability_index(),
                "isoelectric_point": analysed_seq.isoelectric_point(),
                "aliphatic_index": calculate_aliphatic_index(str_protein),
                
                # Skład aminokwasów
                "frac_A": analysed_seq.get_amino_acids_percent().get('A', 0),
                "frac_R": analysed_seq.get_amino_acids_percent().get('R', 0),
                "frac_C": analysed_seq.get_amino_acids_percent().get('C', 0),
                "frac_D": analysed_seq.get_amino_acids_percent().get('D', 0),
                "frac_G": analysed_seq.get_amino_acids_percent().get('G', 0),
                "frac_L": analysed_seq.get_amino_acids_percent().get('L', 0),
                
                "target": 0  # 0 oznacza zbiór negatywny
            }
            results.append(features)
            
        except Exception:
            rejected_err += 1
            continue

    print(f"Zakończono analizę.")
    print(f" - Odrzucono (długość): {rejected_len}")
    print(f" - Odrzucono (błędy biochemiczne): {rejected_err}")
    
    return pd.DataFrame(results)

# --- URUCHOMIENIE ---
input_fasta = "large_negative_dataset.fasta"
output_csv = "ecoli_features_negative_with_seq.csv"

if __name__ == "__main__":
    if not os.path.exists(input_fasta):
        print(f"Błąd: Nie znaleziono pliku {input_fasta}. Wygeneruj go najpierw skryptem 'generate_negatives.py'.")
    else:
        df_negative = calculate_physicochemical_features(input_fasta)
        
        # Zapis do CSV
        df_negative.to_csv(output_csv, index=False)
        print(f"Sukces! Zapisano {len(df_negative)} rekordów do pliku: {output_csv}")
        
        # Podgląd
        cols_to_show = ['id', 'orf_gc', 'orf_entropy', 'molecular_weight']
        print("\nPodgląd nowych cech:")
        print(df_negative[cols_to_show].head())