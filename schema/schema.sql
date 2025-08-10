-- Gene database schema
-- This schema defines the structure for storing gene information and protein atlas expression data
DROP table if exists gene;
CREATE TABLE gene (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    hgnc_symbol TEXT ,
    ensembl_gene_id TEXT not null,
    chromosome_name TEXT,
    description TEXT,
    strand TEXT,
    transcript_count INTEGER,
    gene_start INTEGER,
    gene_end INTEGER,
    unique_exon_count INTEGER,
    gene_length INTEGER,
    UNIQUE(ensembl_gene_id)
);
DROP table if exists protein_atlas_expression;

CREATE TABLE protein_atlas_expression (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_id INTEGER NOT NULL,
    tissue TEXT NOT NULL,
    nTPM REAL,
    FOREIGN KEY(gene_id) REFERENCES gene(id) ON DELETE CASCADE,
    UNIQUE(gene_id, tissue)
);

-- CDoseMap CNV data
DROP TABLE IF EXISTS CDoseMap;
CREATE TABLE CDoseMap (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_id INTEGER NOT NULL,
    cnv_type TEXT,
    cytoband TEXT,
    size INTEGER,
    discovery_sig TEXT,
    known_gd TEXT,
    gnomad_constrained_genes TEXT,
    db_source TEXT,
    article TEXT,
    FOREIGN KEY(gene_id) REFERENCES gene(id) ON DELETE CASCADE
);

-- GTEx gene statistics enriched data
DROP TABLE IF EXISTS gtex_gene_stats;
CREATE TABLE gtex_gene_stats (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_id INTEGER NOT NULL,
    mean REAL,
    median REAL,
    std REAL,
    tau REAL,
    n_zero INTEGER,
    n_below_low INTEGER,
    n_above_low INTEGER,
    n_above_high INTEGER,
    FOREIGN KEY(gene_id) REFERENCES gene(id) ON DELETE CASCADE,
    UNIQUE(gene_id)
);


DROP TABLE IF EXISTS gtex_tissue_tpm;
CREATE TABLE gtex_tissue_tpm (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_id INTEGER NOT NULL,
    tissue TEXT NOT NULL,
    TPM REAL,
    FOREIGN KEY(gene_id) REFERENCES gene(id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS drugs_dgi;
CREATE TABLE drugs_dgi (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_id INTEGER NOT NULL,
    n_drugs INTEGER,
    n_interaction_types INTEGER,
    mean_interaction_score REAL,
    max_interaction_score REAL,
    sum_interaction_score REAL,    
    n_approved_drugs INTEGER,
    has_immunotherapy BOOLEAN,
    has_anti_neoplastic BOOLEAN,
    n_types INTEGER,
    FOREIGN KEY(gene_id) REFERENCES gene(id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS drugs_central;
CREATE TABLE drugs_central (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_id INTEGER NOT NULL,
    n_drugs INTEGER,
    mean_act_value REAL,
    median_act_value REAL,
    n_act_type INTEGER,
    n_target_class INTEGER,
    best_drug_name TEXT,
    best_act_value REAL,
    best_act_type TEXT,
    best_target_class TEXT,
    FOREIGN KEY(gene_id) REFERENCES gene(id) ON DELETE CASCADE,
    UNIQUE(gene_id)
);

-- gnomAD constraint metrics (combined versions v2.1.1 & v4.1)
DROP TABLE IF EXISTS gnomad;
CREATE TABLE gnomad (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_id INTEGER NOT NULL,
    mis_oe_v2_1_1 REAL,
    mis_mu_v2_1_1 REAL,
    mis_pphen_oe_v2_1_1 REAL,
    lof_mu_v2_1_1 REAL,
    lof_pLI_v2_1_1 REAL,
    lof_oe_v2_1_1 REAL,
    mis_z_score_v2_1_1 REAL,
    lof_z_score_v2_1_1 REAL,
    lof_oe_v4_1 REAL,
    lof_pLI_v4_1 REAL,
    lof_z_raw_v4_1 REAL,
    mis_oe_v4_1 REAL,
    mis_z_raw_v4_1 REAL,
    mis_pphen_oe_v4_1 REAL,
    FOREIGN KEY(gene_id) REFERENCES gene(id) ON DELETE CASCADE,
    UNIQUE(gene_id)
);

-- Indexes for better query performance
CREATE INDEX idx_gene_hgnc_symbol ON gene(hgnc_symbol);
CREATE INDEX idx_gene_ensembl_gene_id ON gene(ensembl_gene_id);
CREATE INDEX idx_expr_gene ON protein_atlas_expression(gene_id);
CREATE INDEX idx_expr_tissue ON protein_atlas_expression(tissue);
CREATE INDEX idx_cdose_gene ON CDoseMap(gene_id);
CREATE INDEX idx_cdose_cnv_type ON CDoseMap(cnv_type);
CREATE INDEX idx_cdose_cytoband ON CDoseMap(cytoband);
CREATE INDEX idx_drugs_dgi ON drugs_dgi(gene_id);
CREATE INDEX idx_drugs_central ON drugs_central(gene_id);
CREATE INDEX idx_gnomad_gene ON gnomad(gene_id);
