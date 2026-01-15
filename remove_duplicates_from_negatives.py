import pandas as pd
import numpy as np

# ======================================================================
# KROK 1: PRZYGOTOWANIE I WERYFIKACJA DANYCH
# ======================================================================
# Założenie: 'neg_df' i 'pos_filtered' są już wczytane w poprzednich komórkach.
# Jeśli uruchamiasz to w nowym notebooku, odkomentuj poniższe linie:
# neg_df = pd.read_csv("negatives_diverse_100k.csv")
# pos_filtered = pd.read_csv("positives_with_upstream.csv") # przykładowo


# Sprawdzenie obecności zmiennych
# if 'neg_df' not in locals():
#     raise ValueError("BŁĄD: Zmienna 'neg_df' nie istnieje! Wczytaj najpierw dane negatywne.")
# if 'pos_filtered' not in locals():
#     print("⚠️ Ostrzeżenie: Brak 'pos_filtered'")
# else:
#     target_n = len(pos_filtered)

# print(f"✓ Wczytano zmienną 'neg_df'. Rozmiar oryginalny: {len(neg_df):,} próbek")
# print(f"✓ Cel balansowania (liczba unikalnych pozytywów): {target_n}")

# Sprawdzenie typu danych
# print(f"Typ danych kolumny 'sorf': {neg_df['sorf'].dtype}")
# expected_cols = ["sorf", "upstream_30bp", "accession"]
# missing_cols = [c for c in expected_cols if c not in neg_df.columns]
# if missing_cols:
#     print(f"⚠️ Uwaga: Brakuje kolumn: {missing_cols}")
# else:
#     print(f"✓ Kolumny kluczowe obecne: {expected_cols}")


# ======================================================================
# KROK 2: ANALIZA DUPLIKATÓW (PRZED USUNIĘCIEM)
# ======================================================================
def duplicates_analysis(neg_df):
    # 2.1 Obliczenia
    total_raw = len(neg_df)
    duplicates_count = total_raw - len(neg_df.drop_duplicates(subset=["sorf"]))
    duplicate_pct = (duplicates_count / total_raw) * 100

    print(f"Liczba duplikatów (dokładnych powtórzeń sekwencji): {duplicates_count:,}")
    print(f"Procent duplikatów: {duplicate_pct:.2f}%")

    # 2.2 Top 10
    print("\nTop 10 najczęściej powtarzanych sekwencji sORF:")
    top_dupes = neg_df["sorf"].value_counts().head(10)

    for i, (seq, count) in enumerate(top_dupes.items(), 1):
        # Skracanie sekwencji dla czytelności
        seq_short = seq[:10] + "..." + seq[-3:] if len(str(seq)) > 15 else str(seq)
        print(f"   {i}. Sekwencja #{seq_short:<18} pojawia się {count:>5} razy")

    # 2.3 Statystyki długości
    lengths = neg_df["sorf"].astype(str).str.len()
    print("\nStatystyki długości sekwencji:")
    print(
        f"   Min: {lengths.min()} bp, Max: {lengths.max()} bp, Średnia: {lengths.mean():.1f} bp"
    )


# ======================================================================
# KROK 3: USUWANIE DUPLIKATÓW
# ======================================================================
def duplicates_removal(neg_df):

    # Usuwanie duplikatów
    neg_df_unique = neg_df.drop_duplicates(subset=["sorf"], keep="first").copy()

    # Asercja weryfikująca
    assert (
        len(neg_df_unique) == neg_df["sorf"].nunique()
    ), "BŁĄD KRYTYCZNY: Duplikaty nadal istnieją!"

    print(f"✓ Usunięto duplikaty. Nowy rozmiar zbioru: {len(neg_df_unique):,}")
    print(
        f"✓ Oszczędzono {len(neg_df_unique) / len(neg_df) * 100:.1f}% oryginalnych danych (unikalne reprezentanty)."
    )
    return neg_df_unique


# ======================================================================
# KROK 4: BALANSOWANIE KLAS
# ======================================================================
def balance_negatives(neg_df_unique, target_n):
    if len(neg_df_unique) >= target_n:
        # Mamy wystarczająco unikalnych negatywów
        neg_balanced = neg_df_unique.sample(
            n=target_n, random_state=42, replace=False
        ).copy()
        print(f"✓ Pomyślnie wylosowano {len(neg_balanced)} unikalnych negatywów.")
        print("   (Stosunek 1:1 do pozytywów)")
    else:
        # Mamy za mało unikalnych negatywów
        print(
            f"⚠️ OSTRZEŻENIE: Liczba unikalnych negatywów ({len(neg_df_unique)}) jest mniejsza niż pozytywów ({target_n})!"
        )
        print(
            "   Używam wszystkich dostępnych unikalnych negatywów (zbiór będzie lekko niezbalansowany)."
        )
        neg_balanced = neg_df_unique.copy()
    return neg_balanced


# ======================================================================
# KROK 5: WERYFIKACJA REZULTATÓW
# ======================================================================
def verify_results(neg_balanced, neg_df):

    checks = {
        "Unikalność sekwencji": len(neg_balanced) == neg_balanced["sorf"].nunique(),
        "Brak NaN (sorf)": neg_balanced["sorf"].isnull().sum() == 0,
        "Zgodność kolumn": set(neg_balanced.columns) == set(neg_df.columns),
    }

    all_passed = True
    for name, result in checks.items():
        status = "OK" if result else "BŁĄD"
        print(f"[{status}] {name}")
        if not result:
            all_passed = False

    # Naprawa ewentualnych braków w upstream
    if "upstream_30bp" in neg_balanced.columns:
        nan_upstream = neg_balanced["upstream_30bp"].isnull().sum()
        if nan_upstream > 0:
            print(f"   -> Naprawiono {nan_upstream} pustych wartości w 'upstream_30bp'")
            neg_balanced["upstream_30bp"] = neg_balanced["upstream_30bp"].fillna("")

    if all_passed:
        print("\n✓ SUKCES: Zbiór negatywny 'neg_balanced' jest gotowy do treningu.")
    else:
        print("\n❌ BŁĄD: Znaleziono problemy w danych!")


# ======================================================================
# KROK 6: ZAPIS WYNIKÓW (OPCJONALNIE)
# ======================================================================
# neg_balanced.to_csv("negatives_balanced.csv", index=False)
# neg_df_unique.to_csv("negatives_unique_full.csv", index=False)
