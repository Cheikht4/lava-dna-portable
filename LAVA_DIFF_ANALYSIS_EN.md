# 🧬 Exhaustive Analysis of LAVA-DNA Modifications
## Comparison between `pseudogene/lava-dna` (original) and the modified fork (2026)

**Date of Analysis**: April 26, 2026  
**Original repository**: https://github.com/pseudogene/lava-dna.git  
**Analysis author**: Automated bioinformatics analysis  

---

## 📊 Global Statistical Summary

| Category | Original | Modified Fork | Delta |
|----------|----------|--------------|-------|
| **Main Perl scripts** | 1 (`lava.pl`, 2818 lines) | 2 (`lava_stem_primer.pl` 3934L, `lava_loop_primer.pl` 3452L) | +1 file, +4568 lines |
| **Perl Modules (lib/)** | 14 files | 21 files (+7 new) | +7 modules |
| **Web application (Python)** | 0 | 2 files (1741 lines) | +1741 lines |
| **HTML templates** | 0 | 5 files (1542 lines) | +1542 lines |
| **Deployment** | 0 | 4 files (394 lines) | +394 lines |
| **Documentation** | 1 (`README`) | 6 files | +5 documents |
| **Utility script (`slava.pl`)** | 1 (9808 bytes) | 0 (currently unsupported) | -1 file |

---

## 1. 🔴 DELETED FILES

### 1.1 `lava.pl` → Replaced by two specialized scripts
- **Original role**: Single monolithic LAMP design script (2818 lines)
- **Reason for removal**: Architectural separation into two distinct pipelines (STEM and LOOP) to enable primer design adapted to each LAMP topology
- **Biological impact**: Allows independent processing of STEM primers (positioned between F1c and B1c) and LOOP primers (positioned between F1c/F2 and B1c/B2), with geometric constraints specific to each architecture

### 1.2 `slava.pl` → Removed (Unsupported in V3)
- **Original role**: "Sliding LAVA" script. This wrapper sliced very long genomic sequences (e.g., full 30kb genomes) into overlapping segments (`slava_segment_max_length`) to parallelize or speed up Primer3 analysis.
- **Reason**: Has not been ported to the new Python/Flask architecture yet. Currently, the web interface submits the entire genome at once. (This is a feature that could be reintegrated into the Python backend in the future if full-genome analysis becomes too slow).

---

## 2. 🟢 ADDED FILES — Main Perl Scripts

### 2.1 `lava_stem_primer.pl` [NEW] — 3934 lines (153 KB)
**Nature**: Specialized version derived from `lava.pl` for the STEM architecture

**LAMP STEM Architecture**:
```
5' ←── LAMP Signature (Dynamic Length 'L') ──→ 3'
F3 ──(≥1)── F2 ──(≥1)── F1──(≥1)── FL════BL ──(≥1)── B1 ──(≥1)── B2 ──(≥1)── B3
                         └── Inner+STEM Pair (~40% of L) ──┘
             └──────── Middle Pair (~76% of L) ────────┘
└──────────────────── Outer Pair (100% of L) ──────────────────────┘
```

**Major modifications compared to `lava.pl`**:

| Feature | `lava.pl` (original) | `lava_stem_primer.pl` (fork) |
|---|---|---|
| Mismatch tolerance | None (100% conservation) | `getOligosWithMismatchTolerance()` — position-by-position IUPAC analysis |
| Degeneracy management | Non-existent | Parameters `max_total_degenerate_bases`, `max_consecutive_degenerate_bases`, `max_3prime_degenerate_bases` |
| Signature validation | None | `calculateSignatureIntersection()` — coverage intersection across all primers |
| Results export | Simple text file | Per-signature files with validation report (VALID/REJECT) + amplified/excluded FASTA |
| Combinatorial analysis | None | `analyzeSignatureCombinations()` — search for optimal signature combination for maximum coverage |
| STEM parameters | Non-existent | `stem_primer_target_length/min/max`, `stem_primer_target_tm/min/max` |
| Imported modules | No specific modules | `LLNL::LAVA::Core`, `LLNL::LAVA::Validator` |
| Entropy filtering | Basic | Configurable `entropy_threshold` parameter |

**New functions added**:
- `getOligosWithMismatchTolerance()` — Mismatch tolerance pipeline with IUPAC codes
- `calculateSignatureIntersection()` — Coverage intersection calculation across all primers (6 or 8)
- `analyzeSignatureCombinations()` — Combinatorial analysis to maximize genomic coverage
- `generateCombinations()` — Recursive combination generation
- `createPerSignatureFiles()` — Detailed per-signature export with validation status
- `createAmplificationFiles()` — FASTA export of amplified and excluded sequences

### 2.2 `lava_loop_primer.pl` [NEW] — 3452 lines (140 KB)
**Nature**: Specialized version derived from `lava.pl` for the LOOP architecture

**Key differences compared to `lava_stem_primer.pl`**:

| Aspect | STEM | LOOP |
|--------|------|------|
| Additional primers | FSTEM/BSTEM (between F1c and B1c) | FLOOP/BLOOP (between F1c/F2 and B1c/B2) |
| Coverage parameter | `include_stem_primers` | `include_loop_primers` with threshold `min_signature_coverage` (default 70%) |
| Geometric penalty function | Unified via `PipelineUtils.pm` | Unified via `PipelineUtils.pm` |
| Overlap resolution | `max_overlap_percent`, `resolve_overlap_by` | `max_overlap_percent`, `resolve_overlap_by` (penalty/coverage) |
| Signature validation | Unified via `PipelineUtils.pm` | Unified via `PipelineUtils.pm` |
| Stored signature tags | Unified via `PipelineUtils.pm` | Unified via `PipelineUtils.pm` |

---

## 3. 🟢 ADDED FILES — New Perl Modules

### 3.1 `lib/LLNL/LAVA/Core.pm` [NEW] — 115 lines
**Role**: Shared utility function core between STEM and LOOP

**Exported functions**:
- **`calculate_proportional_geometry($L)`** — Calculates inter-primer target distances based on proportional ratios (F3-F2: 12%, F2-F1: 18%, F1-B1: 40%)
- **`generateSigmoidPenalty($actual, $target, $plateau_ratio, $k_slope)`** — Sigmoid penalty function replacing the original parabola. Comfort zone (plateau 0) at ±25% of target, smooth logistic rise beyond
- **`generateDistancePenalties($maxDistance, $targetLength)`** — Generates a penalty array for all possible distances
- **`countDegenerateBases($sequence)`** — Counts non-standard bases (IUPAC degenerate) in a sequence

**Biological justification**: The sigmoid curve (vs. parabola) more faithfully models the actual thermodynamic behavior of Bst polymerase at 65°C: a wide tolerance zone followed by a smooth increasing penalty, avoiding abrupt rejections of potentially viable candidates.

### 3.2 `lib/LLNL/LAVA/Validator.pm` [NEW] — 449 lines
**Role**: Centralized validation module for LAMP primers and signatures

**Exported functions**:
- **`checkPrimerMismatchTolerance()`** — 4-phase algorithm:
  1. Target region extraction (gap-aware)
  2. SENSE vs ANTISENSE orientation test with automatic reverse-complement
  3. Position-by-position analysis with IUPAC code generation under constraints (max total, max consecutive, max at 3')
  4. Final validation with 3' zone protection and mismatch tolerance
- **`isIUPACCompatible($base, $iupac_code)`** — Verifies IUPAC base-to-code compatibility
- **`rev_comp($seq)`** — Reverse complement with full IUPAC support
- **`generateIUPACCode($bases_ref)`** — Generates the IUPAC code for a set of bases
- **`getPrimerTargetedSequences()`** — Identifies sequences targeted by a given primer
- **`validateCompleteSignatureSpacing()`** — Validates spacing and non-overlap of all primers in a signature

**Biological impact**: This module enables LAMP primer design on highly variable viruses (e.g., Dengue) by tolerating controlled genomic diversity via IUPAC codes, while protecting the 3' end critical for enzymatic extension.

### 3.3 `lib/Lava/Core.pm` [NEW] — 1001 lines
**Role**: Object-oriented refactoring core encapsulating the former `lava.pl` pipeline

**Main functions**:
- `run_lava_loop()` / `run_lava_stem()` — Programmatic entry points (vs. CLI)
- `buildReversePrimers()` — Reverse primer construction via reverse complement
- `analyzeAll()` — Batch analysis via PrimerAnalyzer
- `enumeratePairs()` — Forward/Reverse compatible pair enumeration
- `reducePairInfosByPenalty()` — Best pair selection by score
- `reducePrimersByOverlap()` — Maximum overlap filtering
- `reduceSignaturesByOverlap()` — Complete LAMP signature filtering
- `validateF1B1Spacing()` — F1/B1 spacing validation

### 3.4 `lib/Lava/Enumerator/StemConserved.pm` [NEW] — 81 lines
**Role**: Specialized enumerator for STEM primers, inheriting from `Primer3Conserved`

**Functionality**: Converts named parameters `stem_primer_*` to internal Primer3 tags (`PRIMER_INTERNAL_*`), providing a clear user interface.

### 3.5 `lib/LLNL/LAVA/PipelineUtils.pm` [NEW] — 983 lines
**Role**: Centralization module for STEM and LOOP pipelines (Phase 36 Harmonization)

**Exported Functions**:
- **`calculateSignatureIntersection()`** — Unified logic for signature intersection (6 or 8 primers).
- **`createPerSignatureFiles()`** — Unified export of per-signature results with validation statuses and IUPAC details.
- **`createAmplificationFiles()`** — Unified export of amplified FASTA files.
- **`calculateDynamicPairLengths()`** — Shared dynamic calculation of target inter-primer distances.
- **`reduceSignaturesByOverlap()`** — Harmonization of overlap resolution.
**Impact**: Removed over 1000 lines of duplicated code between STEM and LOOP scripts, ensuring strict algorithmic parity between both approaches.

---

## 4. 🟡 MODIFIED FILES — Existing Perl Modules

### 4.1 `lib/Bio/Tools/Run/Primer3.pm` — 2 modifications

| Modification | Detail |
|---|---|
| Added Primer3 parameters | `PRIMER_TM_FORMULA`, `PRIMER_SALT_CORRECTIONS`, `PRIMER_THERMODYNAMIC_ALIGNMENT`, `SEQUENCE_EXCLUDED_REGION` |
| Result parsing fix | `split '='` → `split('=', $_, 2)` to handle values containing `=` |

**Justification**: The new parameters enable the SantaLucia Tm formula (more accurate for LAMP at 65°C) and SantaLucia salt correction, essential for reliable thermodynamic design.

### 4.2 `lib/LLNL/LAVA/OligoEnumerator/Primer3Conserved.pm` — Major modifications

| Modification | Lines | Description |
|---|---|---|
| **Sliding window entropy** | +50 lines | Window-based smoothing at primer size, marking high-entropy regions as `SEQUENCE_EXCLUDED_REGION` for Primer3 |
| **Zero results diagnostic** | +27 lines | Display of excluded regions and `PRIMER_*_EXPLAIN` when no oligos are found |
| **Strict conservation filter removal** | -25 lines | The former filter requiring 100% identity across all MSA sequences is removed |
| **Post-Primer3 homopolymer filter** | +20 lines | Manual filtering of identical base runs (AAAA, CCCC...) as Primer3 does not always filter them |
| **IUPAC helper** | +17 lines | `_getIUPAC()` function for bases → IUPAC code conversion |
| **IUPAC code disabled (commented)** | +23 lines (commented) | Degenerate sequence generation directly in the enumerator — disabled in favor of the Validator |

**Biological justification**: The removal of the 100% conservation filter is the most important change. The original required each primer to be perfectly identical across ALL MSA sequences — making design on highly variable RNA viruses impossible. The new system delegates diversity management to `Validator.pm` via IUPAC codes.

### 4.3 `lib/LLNL/LAVA/PrimerSet/LAMP.pm` — 1 addition

| Modification | Description |
|---|---|
| **`getStemLocationSummary()`** | +68 lines — New method returning genomic positions of STEM primers (FSTEM/BSTEM) with sense/antisense strand handling |

### 4.4 `lib/LLNL/LAVA/PrimerSetAnalyzer/PCRPair.pm` — 1 modification

| Modification | Description |
|---|---|
| **Explicit import** | +1 line — `use LLNL::LAVA::PrimerSetInfo::PCRPair;` added to ensure module loading |

---

## 5. 🟢 ADDED FILES — Web Interface

### 5.1 `lava_flask_app.py` [NEW] — 1433 lines
**Role**: Complete Flask web application for LAMP design

**Features**:
- FASTA file upload (up to 1 GB)
- Interactive configuration of ~40 LAVA parameters (outer, middle, inner, stem/loop primers)
- Asynchronous Perl script launching via `subprocess`
- Real-time execution monitoring
- Results download (signatures, amplified/excluded FASTA)
- Bilingual interface (French/English) via `TRANSLATIONS` dictionary
- Session authentication
- Past execution history

### 5.2 `launch_lava_smart_kill.py` [NEW] — 308 lines
**Role**: Intelligent launcher with automatic shutdown on browser close (WSL environment)

### 5.3 HTML Templates (5 files, 1542 lines)
| File | Lines | Role |
|---|---|---|
| `templates/base.html` | 128 | Jinja2 base template |
| `templates/index.html` | 592 | Main page — configuration form |
| `templates/monitor.html` | 510 | Real-time execution monitoring |
| `templates/executions.html` | 199 | Execution history |
| `templates/login.html` | 113 | Authentication page |

---

## 6. 🟢 ADDED FILES — Production Deployment

### 6.1 `deployment/` [NEW] — 4 files (394 lines)
| File | Role |
|---|---|
| `deploy.sh` | Automated deployment script (193 lines) |
| `gunicorn_config.py` | Gunicorn configuration (workers, timeouts, bind) |
| `nginx_lava.conf` | Nginx configuration (reverse proxy, SSL, upload limits) |
| `lava-dna.service` | systemd unit for Linux service |

---

## 7. 🟢 ADDED FILES — Documentation

| File | Size | Role |
|---|---|---|
| `README.md` | 10 KB | Main fork documentation (Markdown) |
| `README_Interface.md` | 3.4 KB | Web interface guide |
| `DOCUMENTATION_LAVA.txt` | 13 KB | Detailed technical documentation |
| `LAVA_PARAMETERS_REFERENCE.txt` | 9 KB | Complete reference for ~40 parameters |
| `LAVA_EVOLUTION_JOURNAL.md` | 50 KB | Project evolution journal |
| `Makefile` | 30 KB | Complete build Makefile (vs. `Makefile.PL` only in original) |

---

## 8. 🟡 MODIFIED FILES — Configuration

### 8.1 `.gitignore` — Complete rewrite
- **Original**: Perl build file list (blib/, Makefile, etc.)
- **Fork**: Complete gitignore for Python+Perl project (venv, __pycache__, .DS_Store, logs, IDE, results)

---

## 9. 📋 Summary of Biological and Algorithmic Innovations

### 9.1 Recent Refactoring: Phase 36 Harmonization
The most recent stage of development, "Phase 36," focused on reducing technical debt resulting from the initial split into STEM and LOOP pipelines. By creating the `PipelineUtils.pm` module, we unified the shared logic for signature validation, file I/O, and inter-primer geometry calculation. This refactoring ensures that both STEM and LOOP pipelines maintain perfect algorithmic parity while simplifying maintenance, as updates to validation logic or export formats are now managed in a single codebase.

### 9.2 Genomic Diversity Management
| Aspect | Original | Fork |
|---|---|---|
| Strategy | Strict 100% conservation | IUPAC tolerance with configurable thresholds |
| Target viruses | Conserved only | Highly variable (Dengue, etc.) |
| Genomic coverage | All or nothing | Configurable percentage (default 70%) |
| 3' protection | No | 3' zone protected against mismatches |
| Noise filtering | No | `min_base_frequency` to ignore rare variants |

### 9.3 Thermodynamic Model
| Aspect | Original | Fork |
|---|---|---|
| Spacing penalty | Parabola | Generalized sigmoid (plateau + logistic rise) |
| Tm formula | Basic | SantaLucia (via `PRIMER_TM_FORMULA=1`) |
| Salt correction | Basic | SantaLucia (via `PRIMER_SALT_CORRECTIONS=1`) |
| Geometry | Fixed | Proportional to total signature length |

### 9.4 Software Architecture
| Aspect | Original | Fork |
|---|---|---|
| Main script | 1 monolithic | 2 specialized (STEM + LOOP) |
| Utility modules | 0 | 5 new (`Core.pm`, `Validator.pm`, `Lava::Core`, `StemConserved.pm`, `PipelineUtils.pm`) |
| User interface | CLI only | CLI + Flask Web Interface |
| Deployment | Manual | Automated (Nginx + Gunicorn + systemd) |
| Result Export | Raw text file | Per-signature reports + FASTA + `.params.txt` traceability files |
| Post-design validation | None | Coverage intersection + Combinatorics |

---

## 10. 📁 Comparative File Tree

```
pseudogene/lava-dna (ORIGINAL)          LAVA Fork (MODIFIED)
├── .gitignore                          ├── .gitignore [MODIFIED]
├── .travis.yml                         ├── .travis.yml
├── MANIFEST                            ├── MANIFEST
├── Makefile.PL                         ├── Makefile.PL
├── README                              ├── README
├── environment.yml                     ├── environment.yml
├── lava.pl [REMOVED]                   ├── lava_stem_primer.pl [NEW]
├── slava.pl [REMOVED]                  ├── lava_loop_primer.pl [NEW]
│                                       ├── lava_flask_app.py [NEW]
│                                       ├── launch_lava_smart_kill.py [NEW]
│                                       ├── Makefile [NEW]
│                                       ├── README.md [NEW]
│                                       ├── README_Interface.md [NEW]
│                                       ├── DOCUMENTATION_LAVA.txt [NEW]
│                                       ├── LAVA_PARAMETERS_REFERENCE.txt [NEW]
│                                       ├── LAVA_EVOLUTION_JOURNAL.md [NEW]
│                                       ├── .streamlit/ [NEW]
│                                       ├── deployment/ [NEW]
│                                       │   ├── deploy.sh
│                                       │   ├── gunicorn_config.py
│                                       │   ├── nginx_lava.conf
│                                       │   └── lava-dna.service
│                                       ├── templates/ [NEW]
│                                       │   ├── base.html
│                                       │   ├── index.html
│                                       │   ├── monitor.html
│                                       │   ├── executions.html
│                                       │   └── login.html
│                                       ├── static/ [NEW]
├── lib/                                ├── lib/
│   ├── Bio/Tools/Run/Primer3.pm        │   ├── Bio/Tools/Run/Primer3.pm [MODIFIED]
│   └── LLNL/LAVA/                      │   ├── LLNL/LAVA/
│       ├── Constants.pm                │   │   ├── Constants.pm
│       ├── Oligo.pm                    │   │   ├── Oligo.pm
│       ├── OligoEnumerator.pm          │   │   ├── OligoEnumerator.pm
│       ├── OligoEnumerator/            │   │   ├── OligoEnumerator/
│       │   └── Primer3Conserved.pm     │   │   │   └── Primer3Conserved.pm [MODIFIED]
│       ├── Options.pm                  │   │   ├── Options.pm
│       ├── PrimerAnalyzer.pm           │   │   ├── PrimerAnalyzer.pm
│       ├── PrimerAnalyzer/PCRPrimer.pm │   │   ├── PrimerAnalyzer/PCRPrimer.pm
│       ├── PrimerInfo.pm               │   │   ├── PrimerInfo.pm
│       ├── PrimerSet.pm                │   │   ├── PrimerSet.pm
│       ├── PrimerSet/                  │   │   ├── PrimerSet/
│       │   ├── LAMP.pm                 │   │   │   ├── LAMP.pm [MODIFIED]
│       │   └── PCRPair.pm              │   │   │   └── PCRPair.pm
│       ├── PrimerSetAnalyzer/          │   │   ├── PrimerSetAnalyzer/
│       │   └── PCRPair.pm              │   │   │   └── PCRPair.pm [MODIFIED]
│       ├── PrimerSetInfo.pm            │   │   ├── PrimerSetInfo.pm
│       ├── PrimerSetInfo/PCRPair.pm    │   │   ├── PrimerSetInfo/PCRPair.pm
│                                       │   │   ├── TagHolder.pm
│                                       │   │   ├── Core.pm [NEW]
│                                       │   │   ├── Validator.pm [NEW]
│                                       │   │   ├── PipelineUtils.pm [NEW]
│                                       │   │   └── Lava/ [NEW]
│                                       │   │       ├── Core.pm
│                                       │   │       └── Enumerator/
│                                       │   │           └── StemConserved.pm
├── t/                                  ├── t/
└── t_data/                             └── t_data/
```

---

*Document updated on 04/26/2026 — LAVA-DNA Comparative Bioinformatics Analysis (Including Phase 36 Harmonization)*
