# ============================================================
# config.R — Central configuration for the RNA-seq pipeline
# ============================================================
# Edit this file before running any script.

# ----------------------------
# Output directory
# ----------------------------
BASE_OUT <- "/home/apuhart/Downloads/rna_seq_pip/results_test"

# ----------------------------
# Ensembl mirror
# Options: "useast", "uswest", "asia", "www"
# ----------------------------
ENSEMBL_MIRROR <- "useast"

# ----------------------------
# DESeq2 / filtering
# ----------------------------
MIN_COUNTS      <- 3      # minimum rowSums(counts) to keep a gene
LFC_THRESHOLD   <- 0.3    # |log2FC| threshold for volcano labelling
PADJ_THRESHOLD  <- 0.05   # adjusted p-value threshold

# ----------------------------
# fgsea
# ----------------------------
FGSEA_NPERM     <- 10000
FGSEA_PADJ_CUT  <- 0.3    # pre-filter before dotplot
FGSEA_TOP_N     <- 30     # top N pathways to plot

# MSigDB collections to run  (list of lists: category + optional subcategory)
FGSEA_COLLECTIONS <- list(
  list(category = "H",  subcategory = NULL),
  list(category = "C5", subcategory = "BP"),
  list(category = "C5", subcategory = "CC"),
  list(category = "C5", subcategory = "MF")
)

# ----------------------------
# Heatmap
# ----------------------------
HEATMAP_TOP_N   <- 100    # top N DEGs (by padj) for heatmap

# ----------------------------
# Datasets
# ----------------------------
DATASETS <- list(

  fibro = list(
    samples    = c("C6718", "C6719", "C7522", "F4259", "F4676", "F86"),
    conditions = c("control", "control", "control", "FA", "FA", "FA"),
    sex        = c("M", "F", "F", "M", "M", "F"),
    path       = "/media/apuhart/29d1a4b7-162b-4f05-8e38-da1c77381ffc/SSD_Peio/data_marek/results_nfcore/results_fibro/star_salmon/"
  ),

  iNeurons = list(
    samples    = c("C6718", "C6719", "C7522", "F4259", "F4259ed", "F4676", "F4676ed", "F86ed", "F86_R"),
    conditions = c("control", "control", "control", "FA", "FA_edited", "FA", "FA_edited", "FA_edited", "FA"),
    sex        = c("M", "F", "F", "M", "M", "M", "M", "F", "F"),
    path       = "/media/apuhart/29d1a4b7-162b-4f05-8e38-da1c77381ffc/SSD_Peio/data_marek/results_nfcore/results_iNeurons/star_salmon/"
  ),

  iPSC = list(
    samples    = c("C6718", "C6719", "C7522", "ED86", "ED4259", "ED4676", "F4259", "F4676", "F86"),
    conditions = c("control", "control", "control", "FA_edited", "FA_edited", "FA_edited", "FA", "FA", "FA"),
    sex        = c("M", "F", "F", "F", "M", "M", "M", "M", "F"),
    path       = "/media/apuhart/29d1a4b7-162b-4f05-8e38-da1c77381ffc/SSD_Peio/data_marek/results_nfcore/results_iPSC/star_salmon/"
  )
)

# ----------------------------
# Genes of interest
# ----------------------------
GENES_OF_INTEREST <- unique(c(
  # Fe-S biogenesis — mitochondrial ISC
  "ISCU", "NFS1", "LYRM4", "NDUFAB1", "FXN", "FDX2", "FDXR", "HSPA9", "HSCB", "GLRX5",
  "ABCB7", "ISCA1", "ISCA2", "IBA57", "NFU1",
  # Fe-S biogenesis — CIA pathway
  "CFDP1", "NBP35", "CIAPIN1", "NDOR1", "NARFL", "CIAO1", "CIAO2B", "CIAO2A", "MMS19",
  # Fe-S client proteins — DNA metabolism
  "PRIM2", "DDX11", "DNA2", "BRIP1", "RTEL1", "ERCC2", "POLA1", "POLD1", "POLE",
  "MUTYH", "NTHL1", "KIF4A",
  # Fe-S client proteins — OXPHOS / metabolism
  "NDUFS1", "NDUFS7", "NDUFS8", "NDUFV1", "SDHB", "UQCRFS1", "ACO2", "LIAS", "FECH", "CISD1",
  "ABCE1", "ACO1", "GPAM", "DPYD",
  # ISR / ATF4 axis
  "ATF4", "ATF3", "ASNS", "DDIT3", "GDF15", "TRIB3",
  # One-carbon / serine metabolism
  "PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2", "MTHFD1L", "MTHFD2", "GSS",
  "ALDH1L1", "ALDH1L2",
  # AMPK / mTOR / ferroptosis
  "PRKAA1", "STK11", "MTOR", "GPX4", "SLC25A28", "IREB2", "SLC11A2",
  # NF-kB
  "RELA", "RELB", "NFKB1", "NFKB2", "REL",
  # cGAS-STING
  "TMEM173", "CGAS", "ZBP1", "RIPK1", "RIPK3", "STING1",
  # Cytokines / ECM
  "IL17A", "MMP1", "MMP2", "MMP9",
  # Interferon
  "IRF3", "IRF7", "IRF9", "STAT1", "STAT2", "STAT3",
  # Stress / transcription
  "TP53", "FOXO3", "HIF1A",
  # AP-1
  "JUN", "FOS", "ATF2",
  # Adaptive immunity
  "RORC", "BATF"
))

# ----------------------------
# Neuronal keyword filter (used when neuronal_only = TRUE in fgsea)
# ----------------------------
NEURONAL_KEYWORDS <- c(
  "neuro", "synap", "axon", "dend",
  "glutamate", "gaba", "calcium",
  "neuron", "neuronal", "transmission", "plasticity",
  "mitochond", "oxidative", "respiratory", "electron transport",
  "atp", "energy", "tca", "krebs",
  "iron", "heme", "sulfur", "fe-s", "isc", "cluster",
  "oxidative stress", "ros", "reactive oxygen",
  "glutathione", "redox", "detox",
  "metabolic", "lipid", "fatty acid", "amino acid",
  "one carbon", "folate", "serine",
  "stress response", "unfolded", "upr", "er stress",
  "atf4", "chop", "ddit3", "integrated stress",
  "inflamm", "nfkb", "cytokine", "interferon", "immune",
  "apoptosis", "autophagy", "ferroptosis"
)
