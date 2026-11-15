# Muscle-Contractility-MYOT

# Molecular Mechanisms of Muscle Contractility

* **Student Name:** Jerson Lloyd T. Ortega
* **Course:** Cell and Molecular Biology

## Assignment Information

* **Assignment Number:** 44
* **Gene:** MYOT
* **Official Name:** myotilin
* **Gene Symbol:** MYOT
* **NCBI Gene ID:** 9499
* **Chromosome:** 5q31.2
* **Transcript Accession:** NM_006790.3
* **Variant:** NM_006790.3:c.179C>G
* **Protein Change:** p.Ser60Cys
* **ClinVar Variation ID:** 5836
* **Protein Accession:** NP_006781.1
* **Variant Type:** Missense
* **HGVS Change:**
  * Nucleotide: c.179C>G
  * Protein: p.Ser60Cys
* **Variant Type & Length:** Single nucleotide variant (substitution)
* **Condition:** MYOT-related disorder; Distal myopathy; Myofibrillar myopathy 3
* **Aggregate Classification:** Pathogenic
* **Review Status:** 2 gold stars (criteria provided, multiple submitters, no conflicts)
* **Last Evaluated:** May 7, 2026
* **Date Accessed:** November 15, 2026
* **Submitter:** Clinical laboratories and expert panels via ClinVar
* **Supporting Publications:** PMID: 10942431, 15111675
* **ClinVar Access Note:** Direct web navigation to `/clinvar/variation/5836/` returned a temporary NCBI "Page not found" notice during the NIH website transition on November 15, 2026; live record metadata was verified via the NCBI ClinVar search index (`5836[Variant ID]`) and NCBI E-utilities (`rentrez::entrez_summary(db = "clinvar", id = "5836")`).

## Research Question

How does the missense variant p.Ser60Cys (c.179C>G) in the MYOT gene structurally and functionally affect myotilin's role in sarcomeric Z-disc organization and muscle integrity?

## Methods

* **1. Project Setup:** Install and load required packages (`rentrez` and `ggplot2`), define parameters for Assignment #44 (`MYOT`), and create subfolders (`data`, `results`, `figures`). *(Errors resolved: installed `rentrez` via `install.packages(c("rentrez", "ggplot2"))` in Posit Cloud prior to sourcing, and corrected the `AminoAcid` data-frame column name to match the `ggplot2` aesthetic mapping.)*

```r
install.packages(c("rentrez", "ggplot2"))
library(rentrez)
library(ggplot2)

gene_symbol <- "MYOT"
transcript_accession <- "NM_006790.3"
clinvar_id <- "5836"
variant_protein <- "p.Ser60Cys"

for (x in c("data", "results", "figures")) {
  dir.create(x, showWarnings = FALSE)
}
```

* **2. Gene Retrieval:** Search NCBI Gene for human `MYOT` (`txid9606`) and save summary metadata to `results/gene_summary.txt`.

```r
hits <- entrez_search(
  db = "gene",
  term = paste0(gene_symbol, "[Gene Name] AND txid9606[Organism]"),
  retmax = 20
)
if (length(hits$ids) == 0) stop("Gene not found")

gene_id <- NA_character_
for (id in hits$ids) {
  rec <- entrez_summary(db = "gene", id = id)
  if (toupper(rec$name) == toupper(gene_symbol)) {
    gene_id <- id
    gene_info <- rec
    break
  }
}
if (is.na(gene_id)) stop("Exact human gene not found")

writeLines(c(
  paste("Gene symbol:", gene_info$name),
  paste("NCBI Gene ID:", gene_id),
  paste("Official name:", gene_info$description),
  paste("Chromosome:", gene_info$chromosome)
), "results/gene_summary.txt")
```

* **3. Protein Sequence Retrieval:** Fetch the GenBank record for `NM_006790.3`, extract `/protein_id="NP_006781.1"`, and download the reference FASTA sequence.

```r
gb <- entrez_fetch(
  db = "nuccore", id = transcript_accession,
  rettype = "gb", retmode = "text"
)
pat <- '/protein_id="NP_[0-9]+[.][0-9]+"'
fragment <- regmatches(gb, regexpr(pat, gb))
if (length(fragment) == 0) {
  stop("NP protein ID not found; inspect the transcript record")
}
protein_accession <- gsub('^/protein_id="|"$', "", fragment)
protein_fasta <- entrez_fetch(
  db = "protein", id = protein_accession,
  rettype = "fasta", retmode = "text"
)
if (!startsWith(protein_fasta, ">")) stop("Invalid FASTA")
fasta_file <- file.path("data", paste0(gene_symbol, "_protein.fasta"))
writeLines(protein_fasta, fasta_file)
writeLines(paste("NP accession:", protein_accession),
           "results/protein_accession.txt")
```

* **4. Protein Analysis:** Calculate protein length (498 amino acids), count all 20 standard amino acids, export `results/amino_acid_composition.csv`, and generate the bar chart.

```r
fasta_lines <- readLines(fasta_file)
sequence <- paste(fasta_lines[!grepl("^>", fasta_lines)], collapse = "")
aa <- strsplit(toupper(sequence), "")[[1]]
protein_length <- length(aa)
standard <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
counts <- table(factor(aa, levels = standard))
results <- data.frame(AminoAcid = standard, Count = as.integer(counts))
write.csv(results, "results/amino_acid_composition.csv", row.names = FALSE)

fig1 <- ggplot(results, aes(AminoAcid, Count)) +
  geom_col(fill = "grey45") + theme_bw() +
  labs(title = paste(gene_symbol, "amino-acid composition"),
       x = "Amino acid", y = "Number of residues")
ggsave("figures/01_amino_acid_composition.png", fig1, width = 8, height = 5, dpi = 300)
```

* **5. Variant Review & Biological Visualization:** Query ClinVar Variation ID `5836` and generate the linear protein coordinate map for `p.Ser60Cys`.

```r
if (nzchar(variant_protein)) {
  matched <- regmatches(variant_protein, regexpr("[0-9]+", variant_protein))
  position <- as.integer(matched)
  if (!is.na(position) && position >= 1 && position <= protein_length) {
    fig2 <- ggplot() +
      geom_segment(aes(x = 1, xend = protein_length, y = 0, yend = 0), linewidth = 2) +
      geom_point(aes(x = position, y = 0), size = 4) +
      annotate("text", x = position, y = 0.25, label = variant_protein, size = 4) +
      xlim(0, protein_length) + ylim(-0.6, 0.6) +
      theme_void() +
      labs(title = paste(gene_symbol, "reference protein"),
           x = "Amino-acid position") +
      theme(axis.title.x = element_text())
    ggsave("figures/02_variant_position.png", fig2, width = 8, height = 2.5, dpi = 300)
  }
}
```

## Results and Figures

* **Transcript Accession:** NM_006790.3
* **Retrieved Protein Accession:** NP_006781.1 (Length: 498 amino acids)

## Data Output Files Generated:

* `results/gene_summary.txt`: Contains official gene symbol (`MYOT`), NCBI Gene ID (`9499`), description (`myotilin`), and chromosome location (`5`).
* `results/protein_accession.txt`: Records the confirmed RefSeq protein accession version (`NP_006781.1`).
* `results/amino_acid_composition.csv`: Tabulated counts of all 20 standard amino acids.
* `data/MYOT_protein.fasta`: Reference protein sequence.

## Figures:

1. **Figure 1:** ![Amino Acid Composition Plot](https://drive.google.com/file/d/1ptis24oh193AWD0eC9NTkr8CCiiN5T5V/view?usp=drivesdk)
   * **Interpretation:** Displays counts for all 20 standard amino acids in the MYOT protein, showing that Serine (S, 53 residues), Glutamine (Q, 42 residues), and Alanine (A, 39 residues) are the most frequent residues. This profile describes global amino-acid abundance but does not reveal 3D folding or disease mechanisms.

2. **Figure 2:** ![Variant Position Map](https://drive.google.com/file/d/1-yP0VOM7zJJ7qDgofcRP2lyF9sV0bXyR/view?usp=drivesdk)
   * **Interpretation:** Linear coordinate map of the 498-amino-acid reference protein highlighting the exact location of residue position 60 where the `p.Ser60Cys` missense variant occurs within the N-terminal serine-rich domain.

3. **Figure 3:** ![Protein Structure and Domain Map](https://drive.google.com/file/d/1Ki-S7OIk-NS4d6V-keSzggEFZEEfu3QW/view?usp=drivesdk)
   * **Interpretation:** UniProt feature viewer (`Q9UBF9`) displaying sequence regions, post-translational modification (PTM) sites, the N-terminal serine-rich disordered region (containing Ser60), and the two C-terminal immunoglobulin-like (Ig-like) C2-type domains.

## Variant interpretation (observed effects, clinical association, and predictions)

The *MYOT* `NM_006790.3:c.179C>G` (`p.Ser60Cys`) variant is classified in ClinVar as Pathogenic and is strongly associated with autosomal dominant myofibrillar myopathy 3 (MFM3), distal myopathy, and limb-girdle muscular dystrophy type 1A (LGMD1A). This missense substitution occurs within a well-characterized mutational hotspot in exon 2 and is absent or extremely rare in general population databases.

Clinical and genetic evidence strongly links `p.Ser60Cys` to progressive skeletal muscle disease. Individuals carrying this variant typically present with adult-onset proximal and distal limb weakness, peripheral neuropathy, and cardiomyopathy in a subset of cases. The variant has been repeatedly identified in unrelated families and co-segregates with disease status across affected relatives, supporting a genuine clinical association rather than a benign polymorphism. *MYOT* encodes myotilin, a 498-amino-acid structural protein localized at the sarcomeric Z-disc of striated muscle. Myotilin cross-links filamentous actin (F-actin) into stable bundles and interacts with alpha-actinin, filamin C, and ZASP to anchor thin filaments and transmit contractile force across adjacent sarcomeres. The `p.Ser60Cys` substitution replaces a conserved polar serine residue with a thiol-containing cysteine inside the N-terminal serine-rich domain. Experimental studies provide direct insight into how exon 2 myotilin mutations alter cellular protein homeostasis. Histopathological analyses of skeletal muscle biopsies from patients harboring `p.Ser60Cys` reveal hallmark myofibrillar disorganization, Z-disc streaming, and insoluble sarcoplasmic protein aggregates immunoreactive for myotilin, desmin, and filamin C. Furthermore, biochemical and cell-culture studies of exon 2 hotspot mutants demonstrate that mutant myotilin retains its ability to bind F-actin and alpha-actinin in vitro, but exhibits impaired proteasomal and calpain-mediated degradation, leading to abnormal protein accumulation and defective Z-disc remodeling under mechanical stress. Computational pathogenicity predictors consistently classify `p.Ser60Cys` as deleterious due to evolutionary conservation at position 60, and structural modeling predicts that introducing an unpaired cysteine sulfhydryl group in the disordered N-terminal region may alter local phosphorylation motifs or promote aberrant disulfide cross-linking. However, because the N-terminus lacks a high-resolution experimental crystal structure, the exact atomic-level conformational perturbation caused by `p.Ser60Cys` and the mechanisms underlying variable clinical expressivity remain incompletely resolved. Together, converging clinical, histopathological, biochemical, and computational evidence supports that `p.Ser60Cys` impairs myotilin turnover and destabilizes Z-disc architecture, compromising force transmission during muscle contraction.

## Molecular Causal Chain

Gene (MYOT) $\rightarrow$ Protein (Myotilin) $\rightarrow$ Molecular Function (Z-Disc Actin Cross-Linking & Sarcomere Stabilization) $\rightarrow$ Variant (p.Ser60Cys) $\rightarrow$ Effect (Impaired Protein Turnover, Z-Disc Aggregation, and Loss of Contractile Force Transmission)

## Limitations and Validation Experiment

* **Limitations:** Sequence analysis and linear position mapping only provide primary sequence coordinates; they cannot directly measure in vivo Z-disc mechanics, protein degradation kinetics, or clinical severity variations.
* **Validation Experiment:** A cycloheximide-chase protein turnover assay paired with immunofluorescence confocal microscopy in differentiated C2C12 myotubes expressing wild-type versus `p.Ser60Cys` myotilin to directly quantify mutant protein half-life and Z-disc aggregate accumulation under mechanical stretch.

## Questions to answer

* **1. What does your assigned gene encode, and how does this protein contribute to muscle function?** *MYOT* encodes myotilin, a 498-amino-acid structural Z-disc protein in skeletal and cardiac muscle that bundles F-actin filaments and scaffolds alpha-actinin and filamin C to maintain sarcomere alignment and transmit contractile force.
* **2. What was the exact protein accession and length retrieved with R?** Accession `NP_006781.1` with a total length of `498` amino acids.
* **3. What did the amino-acid composition chart show? What can it NOT tell you?** It showed quantitative distributions of all 20 standard amino acids, with Serine (53), Glutamine (42), and Alanine (39) being most frequent. It cannot determine 3D folding, Z-disc binding interfaces, or the functional consequences of a single missense variant.
* **4. What is your fixed ClinVar variant, and what is its current classification and evidence strength?** Variant `p.Ser60Cys` (`c.179C>G`), ClinVar Variation ID `5836`, classified as **Pathogenic** (last evaluated May 7, 2026) with a **2-star review status** (*criteria provided, multiple submitters, no conflicts*), associated with MYOT-related disorder, distal myopathy, and myofibrillar myopathy 3.
* **5. What is the scientifically supported mechanism linking this variant to muscle function, and what remains uncertain?** Substitution of conserved Ser60 with Cys in the N-terminal serine-rich domain impairs normal myotilin turnover and promotes insoluble myotilin-positive Z-disc aggregates and myofibrillar disorganization. The exact atomic-level structural perturbation in the disordered N-terminus and factors driving variable muscle-group involvement remain uncertain.
* **6. What additional experiment would help test the proposed mechanism?** A cycloheximide-chase immunoblotting assay and immunofluorescence microscopy in mechanically stretched myotubes comparing wild-type and `p.Ser60Cys` myotilin to measure degradation rates and Z-disc aggregation.

## (g) Database Links, Accession Versions, Access Dates, and Citations

* **NCBI Gene:** [MYOT Gene Record (ID: 9499)](https://www.ncbi.nlm.nih.gov/gene/9499) 
* **ClinVar Variant:** [ClinVar Variation 5836](https://www.ncbi.nlm.nih.gov/clinvar/?term=5836%5BVariant+ID%5D) 
* **UniProt Entry:** [Q9UBF9 (MYOT_HUMAN)](https://www.uniprot.org/uniprotkb/Q9UBF9) 
* **Supporting Publications:**
  * Hauser, M. A., Horrigan, S. K., Salmikangas, P., et al. (2000). Myotilin is mutated in limb girdle muscular dystrophy 1A. *Human Molecular Genetics*, 9(14), 2141–2147. [PMID: 10942431](https://pubmed.ncbi.nlm.nih.gov/10942431/)
  * Selcen, D., & Engel, A. G. (2004). Mutations in myotilin cause myofibrillar myopathy. *Neurology*, 62(8), 1363–1371. [PMID: 15111675](https://pubmed.ncbi.nlm.nih.gov/15111675/)
* **Reference Citation:** NCBI RefSeq transcript `NM_006790.3` and protein `NP_006781.1` records.
