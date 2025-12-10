# DNAmusic

[![Open DNAmusic1 in MATLAB Online](https://matlab.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/DNAmusic)
## 🎼 Overview
DNAmusic is a small MATLAB collection for DNA/protein sonification (audification). It includes three routines that implement three published approaches to translate biological sequences into musical structures. Each function retrieves an mRNA record from NCBI (or uses a local GenBank file if available), extracts the CDS and protein translation, and produces audio output.

This repository preserves the original algorithmic intent of each method while providing improved robustness and a progress bar during playback.

## ✨ Included Routines
- DNAmusic1
  Implements the algorithm proposed by Nobuo Munakata. Nucleotide and amino acid mappings are kept fixed according to the original method.
- DNAmusic2
  Implements the approach by Ross D. King and Colin G. Angus (Protein Music).
- DNAmusic3
  Implements the approach by Rie Takahashi and Jeffrey Miller. Codon frequency classes are used to modulate note durations; the codon table is embedded in the code for portability.

## 🧠 Key Ideas
- Audification can reveal patterns that may be less obvious in purely visual inspection.
- The three algorithms provide different musical grammars:
  - DNAmusic1: structured mRNA prelude/main/finale with a combined CDS/protein main theme.
  - DNAmusic2: alternative nucleotide mapping and protein patterns per the cited method.
  - DNAmusic3: protein triads with duration driven by codon frequency classes.

## ⚙️ Requirements
- MATLAB
- Bioinformatics Toolbox (required for genbankread/getgenbank)

## 📦 Installation
1. Clone or download this repository.
2. Add the folder to the MATLAB path.

## 🚀 Usage
Basic usage with defaults:
- DNAmusic1
- DNAmusic2
- DNAmusic3

Specify a different mRNA ID:
- DNAmusic1('NM_XXXXX')
- DNAmusic2('NM_XXXXX')
- DNAmusic3('NM_XXXXX')

Adjust the whole note duration (seconds):
- DNAmusic1('NM_XXXXX', 2)
- DNAmusic2('NM_XXXXX', 1)
- DNAmusic3('NM_XXXXX', 2)

During playback, a progress bar indicates the current section (DNAmusic1/2) or the overall protein triad sequence (DNAmusic3).

## 📝 Notes
- If a file named <mRNAid>.GBK is present in the working directory, the functions will use it before querying NCBI.
- Longer mRNA sequences will require more time and memory to convert into audio.
- The sample rate is set to a historical default of 8192 Hz for consistency with the original implementations.

## 🧾 Citation
If you use this code in academic work, please cite the original algorithms and this repository:
- Munakata N. Algorithm for mapping DNA/protein to music (as implemented in DNAmusic1).
- King RD, Angus CG. Protein Music. CABIOS 1996; 12(3):251-252.
- Takahashi R, Miller J. Method described in Genome Biology 2007; 8(5):405.
- Cardillo G. DNAmusic (MATLAB implementations). GitHub repository.

## 👤 Author
Giuseppe Cardillo  
giuseppe.cardillo.75@gmail.com

## 📄 License
See LICENSE file.
