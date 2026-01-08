"""
Level 1: sORF Detection - Heuristic Baseline
Wymaga: NumPy, Scikit-learn, BioPython
Pamięć: ~200MB | CPU: ~1s na 1000 sekwencji
"""

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_auc_score
from collections import Counter
import re
import pandas as pd

class sORFDetectorLevel1:
    """Heuristic-based sORF detector using sequence features"""
    
    def __init__(self):
        self.model = RandomForestClassifier(n_estimators=100, max_depth=10, random_state=42)
        self.scaler = StandardScaler()
        self.kozak_pattern = re.compile(r'[AG]CCAUGG', re.IGNORECASE)
        
    def extract_gc_content(self, seq):
        """Calculate GC content percentage"""
        if len(seq) == 0:
            return 0
        gc = seq.count('G') + seq.count('C') + seq.count('g') + seq.count('c')
        return gc / len(seq)
    
    def kozak_score(self, upstream_seq):
        """Score Kozak consensus sequence (Purine -3, G at -1)"""
        if len(upstream_seq) < 10:
            return 0
        score = 0
        # Position -3 from AUG: purine (A/G)
        if len(upstream_seq) >= 3 and upstream_seq[-3].upper() in 'AG':
            score += 2
        # Position -1 from AUG: G
        if upstream_seq[-1].upper() == 'G':
            score += 3
        return score / 5  # Normalize to 0-1
    
    def calculate_hexamer_entropy(self, seq):
        """Shannon entropy of hexamer frequency"""
        if len(seq) < 6:
            return 0
        hexamers = [seq[i:i+6] for i in range(len(seq)-5)]
        counter = Counter(hexamers)
        freqs = np.array(list(counter.values())) / len(hexamers)
        entropy = -np.sum(freqs * np.log2(freqs + 1e-10))
        return entropy / np.log2(len(counter) + 1)  # Normalize
    
    def codon_usage_bias(self, seq):
        """Calculate codon usage bias (AT vs GC codons)"""
        if len(seq) < 3:
            return 0
        codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
        at_codons = sum(1 for c in codons if c.count('A') + c.count('T') + 
                       c.count('a') + c.count('t') >= 2)
        return at_codons / len(codons) if codons else 0
    
    def extract_features(self, orf_seq, upstream_seq, downstream_seq):
        """
        Extract features from ORF and flanking regions
        
        Args:
            orf_seq: Open reading frame sequence
            upstream_seq: ~30bp upstream of start codon
            downstream_seq: ~30bp downstream of stop codon
        
        Returns:
            Feature vector
        """
        features = {}
        
        # ORF properties
        features['orf_length'] = len(orf_seq) / 100  # Normalize (typical: 30-150 aa)
        features['orf_gc'] = self.extract_gc_content(orf_seq)
        features['orf_entropy'] = self.calculate_hexamer_entropy(orf_seq)
        features['orf_codon_bias'] = self.codon_usage_bias(orf_seq)
        
        # Upstream context
        features['upstream_gc'] = self.extract_gc_content(upstream_seq)
        features['kozak_score'] = self.kozak_score(upstream_seq)
        features['upstream_codon_bias'] = self.codon_usage_bias(upstream_seq)
        
        # Downstream context
        features['downstream_gc'] = self.extract_gc_content(downstream_seq)
        features['downstream_entropy'] = self.calculate_hexamer_entropy(downstream_seq)
        
        return np.array(list(features.values()))
    
    def fit(self, features, labels):
        """Train the model"""
        self.scaler.fit(features)
        features_scaled = self.scaler.transform(features)
        self.model.fit(features_scaled, labels)
        
    def predict(self, features):
        """Predict sORF probability"""
        features_scaled = self.scaler.transform(features.reshape(1, -1))
        return self.model.predict_proba(features_scaled)[0, 1]
    
    def batch_predict(self, features_list):
        """Predict for multiple sequences"""
        features_scaled = self.scaler.transform(features_list)
        return self.model.predict_proba(features_scaled)[:, 1]


# ============================================================================
# Example usage with synthetic training data
# ============================================================================

if __name__ == "__main__":
    # Generate synthetic training data (replace with real data)
    # np.random.seed(42)
    # n_samples = 500
    
    # smORF (positive): ~50-150 bp
    smorf_features = pd.read_csv("../data/ecoli_only_100000.csv").sample(n=16680, random_state=42, ignore_index=True)
    # smorf_features = np.random.randn(17000 // 2, 9) * 0.3 + \
    #                  np.array([0.5, 0.45, 0.6, 0.5, 0.4, 0.7, 0.45, 0.5, 0.55])
    smorf_labels = np.ones(smorf_features.shape[0])
    
    # Non-smORF (negative): ~200-3000 bp
    non_sorf = pd.read_csv("ecoli_features_negative_with_seq.csv")
    # non_smorf_features = np.random.randn(n_samples // 2, 9) * 0.3 + \
    #                      np.array([2.0, 0.50, 0.3, 0.55, 0.5, 0.3, 0.50, 0.55, 0.40])
    non_smorf_labels = np.zeros(non_sorf.shape[0])
    
    # X = np.vstack([smorf_features['sorf'], non_sorf['dna_sequence']])
    X = pd.concat([smorf_features['sorf'], non_sorf['dna_sequence']], axis=0)
    y = np.hstack([smorf_labels, non_smorf_labels])
    print(f"Training data shape: {X.shape}, Labels shape: {y.shape}")
    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Train model
    detector = sORFDetectorLevel1()
    detector.fit(X_train, y_train)
    
    # Evaluate
    y_pred_proba = detector.batch_predict(X_test)
    y_pred = (y_pred_proba > 0.5).astype(int)
    
    print("="*60)
    print("LEVEL 1: Heuristic sORF Detector - Performance")
    print("="*60)
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred, target_names=['Non-smORF', 'smORF']))
    print(f"ROC-AUC Score: {roc_auc_score(y_test, y_pred_proba):.4f}")
    print("\nFeature Importance:")
    feature_names = ['orf_length', 'orf_gc', 'orf_entropy', 'orf_codon_bias',
                     'upstream_gc', 'kozak_score', 'upstream_codon_bias',
                     'downstream_gc', 'downstream_entropy']
    for name, imp in sorted(zip(feature_names, detector.model.feature_importances_),
                           key=lambda x: x[1], reverse=True):
        print(f"  {name}: {imp:.4f}")
    
    # Example: predict single sequence
    test_features = np.random.randn(1, 9)
    test_features = (test_features - np.mean(X_train, axis=0)) / np.std(X_train, axis=0)
    prob = detector.predict(test_features[0])
    print(f"\nExample prediction: {prob:.4f} (smORF confidence)")
