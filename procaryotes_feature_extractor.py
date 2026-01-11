class EcoliFeatureExtractor:
    """E. coli specific feature extraction"""

    SHINE_DALGARNO = ["AGGAGG", "GGAGG", "GAGG", "GGAG"]

    ECOLI_START_CODONS = {
        "ATG": 1.0,
        "GTG": 0.52,  # Very common in E. coli!
        "TTG": 0.28,
        "CTG": 0.15,
    }

    def extract_features(self, seq, upstream_30bp=""):
        """Extract E. coli features"""
        features = {}

        # 1. SHINE-DALGARNO (BAKTERIA ONLY!)
        rbs_score = 0
        best_distance = 999
        for rbs in self.SHINE_DALGARNO:
            pos = upstream_30bp.rfind(rbs)
            if pos != -1:
                distance = 30 - (pos + len(rbs))
                if 5 <= distance <= 12:
                    rbs_score = max(rbs_score, (13 - distance) / 8.0)
                    best_distance = min(best_distance, distance)

        features["rbs_score"] = rbs_score
        features["rbs_distance"] = best_distance / 20.0  # Normalize

        # 2. E. coli start codon preference
        start_codon = seq[:3].upper()
        features["start_codon_strength"] = self.ECOLI_START_CODONS.get(start_codon, 0)

        # 3. AT bias (E. coli specific)
        at_count = seq.count("A") + seq.count("T")
        features["at_bias"] = at_count / len(seq)

        # 4. ORF length (sORF: 30-150 bp)
        features["orf_length_norm"] = min(len(seq) / 150.0, 1.0)

        # 5. GC content
        gc_count = seq.count("G") + seq.count("C")
        features["gc_content"] = gc_count / len(seq)

        return list(features.values())
