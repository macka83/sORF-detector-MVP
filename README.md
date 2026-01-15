# sORF-detector
**E.coli sORF simplified binary classifier**

---

## 💻 Environment
* **Platform:** WSL Ubuntu
* **Conda Version:** 25.11.1
* **Python Version:** 3.13.11
* **Modules:** Specified in `requirements.txt`

---

## 📊 Input Data

### Positive sORF Dataset
* **Source:** [sORFdb](https://sorfdb.computational.bio/download)
* **File:** `sorfdb.tsv.gz` (E.coli sample)

### Negative sORF Dataset
* **Source:** NCBI
* **Method:** Multiple E.coli strain sequences obtained via programmatic access.

---

## ⚙️ Data Preprocessing
1. **Context Enrichment:** Addition of the upstream sequence (30bp) before the sequence of interest.
2. **Negative Dataset Cleaning:** Removal of exact duplicates from the negative dataset to prevent data leakage.
3. **Positive Dataset Redundancy Reduction:** Removal of duplicated sequences with >90% similarity using CD-HIT.

---

## 🧬 Numeric Feature Extraction
The following features are extracted for the classification task:

* **RBS Score:** Strength of the Shine-Dalgarno sequence.
* **Start Codon Strength:** Weighted importance of the start codon (ATG/GTG/TTG).
* **Sequence Length Normalisation:** Scaling the length of the ORF.
* **GC Content:** Percentage of G and C nucleotides in the ORF.
* **AT Bias:** AT nucleotide bias analysis.
* **GC Upstream:** GC content in the 30bp upstream region.

---

## 🤖 Training & Model
* **Scaling:** `StandardScaler` (Applied to numeric features; validity for specific sORF biological distributions to be verified).
* **Algorithm:** `RandomForestClassifier`

---

## Files Description
* **fetch_upstream_from_ncbi_for_pos_and_neg.ipynb** add upstream sequence to positive sORF  
* **get_ecoli_source.ipynb** get multiple E.coli sequences from NCBI
* **generate_sORF_negative_sequences.py** read multiple E.coli sequences and generate negative sORF dataset 

* **Feature_extraction_and_training.ipynb** apply Numeric Feature Extraction and Training
