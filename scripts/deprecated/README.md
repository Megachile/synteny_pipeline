# Deprecated Scripts

This folder contains old versions of scripts that have been superseded.

## 04_blast_targets_DEPRECATED_20251103.py

**Deprecated:** 2025-11-03
**Reason:** Duplicate implementation. The "optimized" version was the canonical one actually used.
**Replaced by:** `scripts/04_blast_targets.py` (formerly `04_blast_targets_optimized.py`)

**Key differences:**
- Old version had more complex Exonerate extraction logic embedded
- New version is cleaner, more modular
- Both had the same clustering bug (not respecting reading frames) which was fixed in canonical version

**Action taken:** Renamed optimized version to canonical name, deprecated this duplicate.
