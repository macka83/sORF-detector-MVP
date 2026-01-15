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

# Proporcje zbioru negatywnego (wg analizy)
RATIOS = {
    "intergenic": 0.60,  # Naturalne tło
    "frameshift": 0.20,  # "Trudne" negatywy (przesunięta ramka)
    "shuffled": 0.20,  # Szum statystyczny
}

# Parametry sORF
MIN_BP = 30  # 10 aa
MAX_BP = 300  # 100 aa (dłuższe traktujemy jako geny referencyjne)
UPSTREAM_LEN = 30

# Stałe taksonomiczne (aby pasowały do positive dataset)
TAXONOMY_TEMPLATE = {
    "phylum": "Pseudomonadota",
    "class": "Gammaproteobacteria",
    "order": "Enterobacterales",
    "family": "Enterobacteriaceae",
    "genus": "Escherichia",
    "species": "Escherichia coli",
}


def get_physiochem_properties(aa_seq):
    """Oblicza parametry białka (kopia logiki z positives)"""
    feats = {
        "molecular-weight": 0,
        "gravy": 0,
        "aromaticity": 0,
        "instability": 0,
        "isoelectric-point": 0,
        "aliphatic-index": 0,
        "boman": 0,
    }
    clean = aa_seq.replace("*", "").replace("X", "")
    if not clean:
        return feats
    try:
        a = ProtParam.ProteinAnalysis(clean)
        feats["molecular-weight"] = a.molecular_weight()
        feats["gravy"] = a.gravy()
        feats["aromaticity"] = a.aromaticity()
        feats["instability"] = a.instability_index()
        feats["isoelectric-point"] = a.isoelectric_point()

        # Prosty Aliphatic Index
        A, V, I, L = (
            clean.count("A"),
            clean.count("V"),
            clean.count("I"),
            clean.count("L"),
        )
        feats["aliphatic-index"] = (A + 2.9 * V + 3.9 * (I + L)) / len(clean) * 100
    except:
        pass
    return feats


def extract_upstream(full_seq, start_idx, strand):
    """Pobiera prawdziwy upstream z genomu"""
    if strand == 1:
        if start_idx >= UPSTREAM_LEN:
            return str(full_seq[start_idx - UPSTREAM_LEN : start_idx])
        return str(full_seq[:start_idx]).rjust(UPSTREAM_LEN, "N")
    else:  # strand -1
        # Dla nici minus, upstream jest "za" ORF na nici plus
        # ORF: [stop ... start] <--- kierunek transkrypcji
        # Genom: 5' ... [RC_ORF] [Upstream_RC] ... 3'
        # Logika: bierzemy fragment ZA (start_idx to początek ORF na nici +)
        # Uwaga: start_idx tutaj to początek matchu w stringu, czyli koniec ORF dla nici -
        # Upraszczając dla Seq object:
        # Pobieramy fragment z genomu i robimy RC
        end_idx = start_idx + UPSTREAM_LEN  # To tylko przybliżenie dla uproszczenia
        # Lepiej: W funkcji głównej operujemy na gotowych stringach
        return "N" * UPSTREAM_LEN  # Placeholder, logika w pętli głównej jest pewniejsza


def generate_entry(seq_id, seq_dna, upstream, source_type, strain_name, assembly):
    """Tworzy wiersz słownika pasujący do formatu CSV"""

    # Translacja
    seq_obj = Seq(seq_dna)
    # Dla sORF zakładamy, że podajemy sekwencję kodującą 5'->3'
    protein = str(seq_obj.translate())

    feats = get_physiochem_properties(protein)

    entry = TAXONOMY_TEMPLATE.copy()
    entry.update(
        {
            "id": seq_id,
            "source": f"synthetic_{source_type}",
            "assembly": assembly,
            "accession": assembly,  # Używamy nazwy pliku jako accession
            "strain": strain_name,
            "sorf": seq_dna,
            "slen": len(seq_dna),
            "start-codon": seq_dna[:3],
            "protein": protein,
            "plen": len(protein),
            "product": "hypothetical protein",
            "rbs": 0,  # Negatyw
            "pfam-hits": "",
            "upstream_30bp": upstream,
            **feats,
        }
    )
    return entry


def main():
    # 1. Znajdź pliki
    if not os.path.exists(INPUT_DIR):
        print("Brak folderu z danymi!")
        return

    files = glob.glob(os.path.join(INPUT_DIR, "*.fasta")) + glob.glob(
        os.path.join(INPUT_DIR, "*.fna")
    )
    num_files = len(files)
    print(f"Znaleziono {num_files} plików genomowych.")

    if num_files == 0:
        return

    # 2. Oblicz ile sekwencji pobrać z każdego pliku (Load Balancing)
    per_file_target = int(np.ceil(TARGET_TOTAL / num_files))
    print(f"Cel: ~{per_file_target} sekwencji na plik (Total: {TARGET_TOTAL})")

    negatives = []

    # Pasek postępu
    pbar = tqdm(total=TARGET_TOTAL, desc="Generowanie zróżnicowanych negatywów")

    for fpath in files:
        if len(negatives) >= TARGET_TOTAL:
            break

        fname = os.path.basename(fpath)
        assembly_id = fname.split(".")[0]

        try:
            # Wczytaj rekord (zakładamy 1 kontig/chromosom dla szybkości)
            record = next(SeqIO.parse(fpath, "fasta"))
            full_seq = str(record.seq).upper()
            seq_len = len(full_seq)

            # Liczniki dla tego pliku
            file_collected = 0

            # --- STRATEGIA 1: INTERGENIC (Pseudo-sORF) ---
            # Szukamy "dziur" w genomie. Losujemy pozycję i szukamy startu.
            # Limit prób, żeby nie utknąć
            for _ in range(50):
                if file_collected >= per_file_target * RATIOS["intergenic"]:
                    break

                # Losowa pozycja
                pos = random.randint(0, seq_len - MAX_BP - 50)
                chunk = full_seq[pos : pos + MAX_BP + 50]

                # Szukamy startu (ATG/GTG/TTG)
                # Uproszczenie: szukamy pierwszego startu w chunku
                start_match = -1
                for start_codon in ["ATG", "GTG", "TTG"]:
                    idx = chunk.find(start_codon)
                    if idx != -1:
                        start_match = idx
                        break

                if start_match != -1:
                    # Szukamy stopu in-frame
                    sub_chunk = chunk[start_match:]
                    for i in range(0, len(sub_chunk) - 3, 3):
                        codon = sub_chunk[i : i + 3]
                        if codon in ["TAA", "TAG", "TGA"]:
                            cand_len = i + 3
                            # Sprawdzamy długość (musi być sORF)
                            if MIN_BP <= cand_len <= 200:  # 200bp max dla sORF
                                # Sprawdźmy czy to nie jest fragment długiego genu (heurystyka)
                                # Jeśli znaleźliśmy go losowo, a nie ma GFF, zakładamy że to intergenic
                                # o ile nie jest bardzo długi.

                                # Pobieramy upstream
                                abs_start = pos + start_match
                                upstream = ""
                                if abs_start >= UPSTREAM_LEN:
                                    upstream = full_seq[
                                        abs_start - UPSTREAM_LEN : abs_start
                                    ]
                                else:
                                    upstream = "N" * UPSTREAM_LEN

                                entry = generate_entry(
                                    f"neg_int_{assembly_id}_{abs_start}",
                                    sub_chunk[:cand_len],
                                    upstream,
                                    "intergenic",
                                    record.description.split(",")[
                                        0
                                    ],  # Próba wyjęcia nazwy szczepu
                                    assembly_id,
                                )
                                negatives.append(entry)
                                file_collected += 1
                                pbar.update(1)
                                break  # Sukces w tej iteracji losowania

            # --- STRATEGIA 2: FRAMESHIFT (Przesunięte ramki) ---
            # Szukamy długiego genu (>500bp) i przesuwamy ramkę
            # Znajdujemy długi ORF
            found_long = False
            for _ in range(20):  # Próby znalezienia długiego ORF
                if file_collected >= per_file_target * (
                    RATIOS["intergenic"] + RATIOS["frameshift"]
                ):
                    break

                pos = random.randint(0, seq_len - 1000)
                chunk = full_seq[pos : pos + 1000]

                # Uproszczone szukanie długiego ORF (start... brak stopu przez długi czas)
                # Tutaj po prostu bierzemy losowy kawałek DNA, który wygląda na kodujący
                # (np. wysokie GC, brak stopów).
                # Lepsza metoda: bierzemy losowy fragment i jeśli nie ma stopów przez 100bp,
                # przesuwamy go o 1 nt i ucinamy sORF.

                # Symulacja Frameshift: Bierzemy fragment, wstawiamy przesunięcie
                # Bierzemy losowy 150bp fragment
                cand_seq = chunk[:150]
                # Sprawdzamy czy ma stopy w ramce 0
                stops = [
                    cand_seq[i : i + 3] in ["TAA", "TAG", "TGA"]
                    for i in range(0, len(cand_seq), 3)
                ]

                if not any(stops):  # Wygląda jak fragment genu
                    # Robimy Frameshift: ucinamy 1 nt z początku
                    shifted_seq = cand_seq[1:]
                    # Musimy znaleźć nowy sztuczny stop lub uciąć
                    # Dla bezpieczeństwa bierzemy po prostu ten fragment jako "coding-like but wrong"

                    # Upstream: Prawdziwy upstream z genomu (przesunięty o 1)
                    abs_start = pos + 1
                    upstream = full_seq[abs_start - UPSTREAM_LEN : abs_start]

                    entry = generate_entry(
                        f"neg_fs_{assembly_id}_{pos}",
                        shifted_seq,
                        upstream,
                        "frameshift",
                        "synthetic",
                        assembly_id,
                    )
                    negatives.append(entry)
                    file_collected += 1
                    pbar.update(1)

            # --- STRATEGIA 3: SHUFFLED (Tasowane) ---
            # Bierzemy to co już zebraliśmy (intergenic) i tasujemy
            needed_shuffled = per_file_target - file_collected
            for _ in range(needed_shuffled):
                if not negatives:
                    break  # Safety

                # Losujemy wzorzec z już zebranych (aby zachować długość i skład)
                template = random.choice(negatives)
                seq_list = list(template["sorf"])
                random.shuffle(seq_list)
                shuffled_seq = "".join(seq_list)

                # Upstream też tasujemy, żeby model nie uczył się korelacji
                up_list = list(template["upstream_30bp"])
                random.shuffle(up_list)
                shuffled_upstream = "".join(up_list)

                entry = generate_entry(
                    f"neg_shuf_{assembly_id}_{random.randint(0,9999)}",
                    shuffled_seq,
                    shuffled_upstream,
                    "shuffled",
                    "synthetic",
                    assembly_id,
                )
                negatives.append(entry)
                pbar.update(1)

        except Exception as e:
            # print(f"Błąd w {fname}: {e}")
            continue  # Skip corrupted files

    pbar.close()

    # 3. Zapis i czyszczenie
    df = pd.DataFrame(negatives)

    # Upewniamy się że kolumny są w dobrej kolejności (jak w positives)
    cols = [
        "id",
        "source",
        "assembly",
        "accession",
        "protein-id",
        "uid",
        "entry-name",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "strain",
        "sorf",
        "slen",
        "start-codon",
        "protein",
        "plen",
        "product",
        "rbs",
        "pfam-hits",
        "gravy",
        "aromaticity",
        "molecular-weight",
        "instability",
        "isoelectric-point",
        "aliphatic-index",
        "boman",
        "upstream_30bp",
    ]

    # Wypełnij brakujące
    for c in cols:
        if c not in df.columns:
            df[c] = ""

    # Przytnij do dokładnej liczby i zapisz
    df = df[cols].head(TARGET_TOTAL)
    df.to_csv(OUTPUT_FILE, index=False)

    print("\n" + "=" * 50)
    print(f"Gotowe! Wygenerowano {len(df)} negatywnych sekwencji.")
    print(f"Plik wynikowy: {OUTPUT_FILE}")
    print(f"Struktura: ~60% Intergenic, ~20% Frameshift, ~20% Shuffled")

    # Szybki check balansu
    print("\nRozkład typów (source):")
    print(df["source"].value_counts())


if __name__ == "__main__":
    main()
