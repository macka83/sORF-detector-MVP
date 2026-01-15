import os
import glob
import random
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import ProtParam
from tqdm import tqdm

# --- KONFIGURACJA ---
INPUT_DIR = "data/ecoli_genomes"
OUTPUT_FILE = "negatives_diverse_13k.csv"
TARGET_TOTAL = 13000

# Proporcje zbioru
RATIOS = {
    "intergenic": 0.60,  # Losowe ORF z tła
    "frameshift": 0.20,  # Przesunięta ramka (trudne negatywy)
    "shuffled": 0.20,  # Całkowity szum (zachowany skład aminokwasowy)
}

# Parametry sORF
MIN_BP = 30  # 10 aa
MAX_BP = 300  # 100 aa
UPSTREAM_LEN = 30


def get_physiochem_properties(aa_seq: str) -> dict:
    """Oblicza parametry fizykochemiczne białka."""
    feats = {
        "molecular_weight": 0.0,
        "gravy": 0.0,
        "aromaticity": 0.0,
        "instability": 0.0,
        "isoelectric_point": 0.0,
        "aliphatic_index": 0.0,
        "boman": 0.0,
    }
    # Usuń stop codon i znaki niejednoznaczne dla obliczeń
    clean = (
        aa_seq.replace("*", "")
        .replace("X", "")
        .replace("J", "")
        .replace("B", "")
        .replace("Z", "")
    )

    if not clean:
        return feats

    try:
        a = ProtParam.ProteinAnalysis(clean)
        feats["molecular_weight"] = a.molecular_weight()
        feats["gravy"] = a.gravy()
        feats["aromaticity"] = a.aromaticity()
        feats["instability"] = a.instability_index()
        feats["isoelectric_point"] = a.isoelectric_point()

        # Aliphatic Index
        A, V, I, L = (
            clean.count("A"),
            clean.count("V"),
            clean.count("I"),
            clean.count("L"),
        )
        if len(clean) > 0:
            feats["aliphatic_index"] = (A + 2.9 * V + 3.9 * (I + L)) / len(clean) * 100
    except Exception:
        # W przypadku błędu (np. bardzo krótka sekwencja lub dziwne znaki) zwracamy zera
        pass

    return feats


def extract_genome_upstream(full_seq: str, start_index: int, length: int = 30) -> str:
    """Pobiera sekwencję upstream z genomu, obsługując krawędzie."""
    if start_index >= length:
        return full_seq[start_index - length : start_index]
    else:
        # Padding 'N' jeśli jesteśmy na początku genomu
        padding = "N" * (length - start_index)
        return padding + full_seq[:start_index]


def create_entry(seq_id, dna_seq, upstream_seq, source_type, assembly_id):
    """Tworzy ujednolicony słownik danych dla pojedynczego rekordu."""

    # Translacja
    protein_seq = str(Seq(dna_seq).translate())

    # Cechy białka
    feats = get_physiochem_properties(protein_seq)

    entry = {
        "id": seq_id,
        "source": f"synthetic_{source_type}",
        "assembly": assembly_id,
        "sorf": dna_seq,
        "slen": len(dna_seq),
        "protein": protein_seq,
        "plen": len(protein_seq),
        "upstream_30bp": upstream_seq,
        # Rozpakowanie cech
        **feats,
    }
    return entry


def generate_intergenic(full_seq, assembly_id, count_needed):
    """Generuje negatywy typu 'intergenic' (przypadkowe ORFy)."""
    generated = []
    seq_len = len(full_seq)
    attempts = 0
    max_attempts = count_needed * 50  # Limit prób, aby nie utknąć

    while len(generated) < count_needed and attempts < max_attempts:
        attempts += 1
        pos = random.randint(0, seq_len - MAX_BP - 50)
        chunk = full_seq[pos : pos + MAX_BP + 50]

        # Szukamy startu
        start_match = -1
        for start_codon in ["ATG", "GTG", "TTG"]:
            idx = chunk.find(start_codon)
            if idx != -1:
                start_match = idx
                break

        if start_match != -1:
            # Szukamy stopu in-frame
            sub_chunk = chunk[start_match:]
            found_stop = False
            valid_len = 0

            for i in range(0, len(sub_chunk) - 3, 3):
                codon = sub_chunk[i : i + 3]
                if codon in ["TAA", "TAG", "TGA"]:
                    valid_len = i + 3
                    found_stop = True
                    break

            if found_stop and (MIN_BP <= valid_len <= MAX_BP):
                abs_start = pos + start_match
                upstream = extract_genome_upstream(full_seq, abs_start, UPSTREAM_LEN)

                entry = create_entry(
                    f"neg_int_{assembly_id}_{abs_start}",
                    sub_chunk[:valid_len],
                    upstream,
                    "intergenic",
                    assembly_id,
                )
                generated.append(entry)

    return generated


def generate_frameshift(full_seq, assembly_id, count_needed):
    """Generuje negatywy przez przesunięcie ramki w regionach kodujących."""
    generated = []
    seq_len = len(full_seq)
    attempts = 0
    max_attempts = count_needed * 50

    while len(generated) < count_needed and attempts < max_attempts:
        attempts += 1
        # Pobieramy losowy fragment, który wygląda na gen (brak stopów przez dłuższą chwilę)
        pos = random.randint(0, seq_len - 1000)
        # Bierzemy fragment i sprawdzamy czy ma stopy
        cand_len = 150  # Długość testowa
        cand_seq = full_seq[pos : pos + cand_len]

        stops = [
            cand_seq[i : i + 3] in ["TAA", "TAG", "TGA"]
            for i in range(0, len(cand_seq), 3)
        ]

        if not any(stops):
            # Wygląda na fragment genu. Przesuwamy o 1bp
            shifted_seq = cand_seq[1:]

            # Pobieramy upstream dla przesuniętego startu
            abs_start = pos + 1
            upstream = extract_genome_upstream(full_seq, abs_start, UPSTREAM_LEN)

            entry = create_entry(
                f"neg_fs_{assembly_id}_{abs_start}",
                shifted_seq,
                upstream,
                "frameshift",
                assembly_id,
            )
            generated.append(entry)

    return generated


def generate_shuffled(existing_negatives, count_needed, assembly_id):
    """Generuje negatywy przez tasowanie sekwencji (zachowuje skład aminokwasowy)."""
    generated = []
    if not existing_negatives:
        return []

    for _ in range(count_needed):
        template = random.choice(existing_negatives)

        # Tasowanie nukleotydów
        seq_list = list(template["sorf"])
        random.shuffle(seq_list)
        shuffled_seq = "".join(seq_list)

        # Tasowanie upstreamu (żeby usunąć sygnał RBS)
        up_list = list(template["upstream_30bp"])
        random.shuffle(up_list)
        shuffled_upstream = "".join(up_list)

        entry = create_entry(
            f"neg_shuf_{assembly_id}_{random.randint(1000,99999)}",
            shuffled_seq,
            shuffled_upstream,
            "shuffled",
            assembly_id,
        )
        generated.append(entry)

    return generated


def main():
    if not os.path.exists(INPUT_DIR):
        print(f"Błąd: Brak folderu {INPUT_DIR}")
        return

    files = glob.glob(os.path.join(INPUT_DIR, "*.fasta")) + glob.glob(
        os.path.join(INPUT_DIR, "*.fna")
    )
    num_files = len(files)
    print(f"Znaleziono {num_files} plików genomowych.")

    if num_files == 0:
        return

    # Load balancing
    per_file_target = int(np.ceil(TARGET_TOTAL / num_files))

    all_negatives = []
    pbar = tqdm(total=TARGET_TOTAL, desc="Generowanie negatywów")

    for fpath in files:
        if len(all_negatives) >= TARGET_TOTAL:
            break

        fname = os.path.basename(fpath)
        assembly_id = fname.split(".")[0]

        try:
            record = next(SeqIO.parse(fpath, "fasta"))
            full_seq = str(record.seq).upper()

            # 1. Intergenic
            n_inter = int(per_file_target * RATIOS["intergenic"])
            inter_data = generate_intergenic(full_seq, assembly_id, n_inter)
            all_negatives.extend(inter_data)
            pbar.update(len(inter_data))

            # 2. Frameshift
            n_frame = int(per_file_target * RATIOS["frameshift"])
            frame_data = generate_frameshift(full_seq, assembly_id, n_frame)
            all_negatives.extend(frame_data)
            pbar.update(len(frame_data))

            # 3. Shuffled (bazuje na już zebranych)
            n_shuf = per_file_target - len(inter_data) - len(frame_data)
            if n_shuf > 0 and (inter_data or frame_data):
                pool = inter_data + frame_data
                shuf_data = generate_shuffled(pool, n_shuf, assembly_id)
                all_negatives.extend(shuf_data)
                pbar.update(len(shuf_data))

        except Exception as e:
            print(f"Pominięto plik {fname}: {e}")
            continue

    pbar.close()

    # Zapis do CSV - tylko istotne kolumny
    # Usunięto: taxonomy, product, uid, etc.
    final_cols = [
        "id",
        "source",
        "assembly",
        "sorf",
        "slen",
        "protein",
        "plen",
        "upstream_30bp",
        "molecular_weight",
        "gravy",
        "aromaticity",
        "instability",
        "isoelectric_point",
        "aliphatic_index",
        "boman",
    ]

    df = pd.DataFrame(all_negatives)

    # Uzupełnienie brakujących kolumn (na wypadek gdyby jakieś obliczenia się nie udały)
    for c in final_cols:
        if c not in df.columns:
            df[c] = None

    # Ostateczny filtr i zapis
    df = df[final_cols].head(TARGET_TOTAL)
    df.to_csv(OUTPUT_FILE, index=False)

    print(f"\nGotowe! Zapisano {len(df)} rekordów do {OUTPUT_FILE}")
    print(df["source"].value_counts())


if __name__ == "__main__":
    main()
