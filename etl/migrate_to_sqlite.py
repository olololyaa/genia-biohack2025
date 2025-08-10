from pathlib import Path
import pandas as pd
import sys
from utils import load_table, create_temp_table_and_insert, check_columns,  cleanup_temp_table, create_database, DB_PATH

# Database and file paths
GENE_CSV = Path("gene_info/genes_info.csv")
CONS_TSV = Path("gene_info/rna_tissue_consensus.tsv")

def migrate_genes(cur):
    """Migrate gene data from CSV to database using temporary table."""
    # TODO refactor, fetch all data from Biomart within a one call
    print("Migrating gene data...")
    
    # Load gene data with validation
    genes_df = load_table(GENE_CSV)
    genes_metrics = load_table(Path("gene_info/biomart_basic_metrics.csv"))
    #drop Strand column from genes_metrics
    genes_metrics = genes_metrics.drop(columns=['Strand'], errors='ignore')
    #import pdb; pdb.set_trace()
    rename_map = {
        'Strand': 'strand',
        'Transcript count': 'transcript_count',
        'Gene start (bp)': 'gene_start',
        'Gene end (bp)': 'gene_end',
        'Chromosome/scaffold name': 'chromosome_name',
        'Unique exon count': 'unique_exon_count',
        'Gene length (bp)': 'gene_length',
        'Gene stable ID': 'ensembl_gene_id'
    }
    #rename gene_metrics columns
    genes_metrics = genes_metrics.rename(columns=rename_map).drop_duplicates()
    #merge 
    genes_all = genes_df.merge(genes_metrics, on='ensembl_gene_id', how='inner')
    #import pdb; pdb.set_trace()
    # find rows where ensembl_gene_id are non unique
    non_unique_genes = genes_all[genes_all.duplicated(subset='ensembl_gene_id', keep=False)]
    #drop one row from genes_all,   keep first occurrence
    genes_all = genes_all[~genes_all.duplicated(subset='ensembl_gene_id', keep='first')]  #bad practice! needs refactoring, 
    #save to csv:
    #non_unique_genes.to_csv("gene_info/non_unique_genes.csv", index=False) # drop duplicates doesn't work
    # Create temporary table
    cur.execute("""
        CREATE TEMPORARY TABLE temp_genes (
            hgnc_symbol TEXT,
            ensembl_gene_id TEXT,
            description TEXT,
            chromosome_name TEXT,
            strand TEXT,
            transcript_count INTEGER,
            gene_start INTEGER,
            gene_end INTEGER,
            unique_exon_count INTEGER,
            gene_length INTEGER
        )
    """)

    

    # Insert all data into temporary table
    temp_data = []
    for row in genes_all.itertuples(index=False):
        temp_data.append((
            getattr(row, 'hgnc_symbol', None),
            getattr(row, 'ensembl_gene_id', None), 
            getattr(row, 'description', None),
            getattr(row, 'chromosome_name', None),
            getattr(row, 'strand', None),
            getattr(row, 'transcript_count', None),
            getattr(row, 'gene_start', None),
            getattr(row, 'gene_end', None),
            getattr(row, 'unique_exon_count', None),
            getattr(row, 'gene_length', None)
        ))

    #import pdb; pdb.set_trace()
    
    cur.executemany(
        """
        INSERT INTO temp_genes (
            hgnc_symbol,
            ensembl_gene_id,
            description,
            chromosome_name,
            strand,
            transcript_count,
            gene_start,
            gene_end,
            unique_exon_count,
            gene_length
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        temp_data
    )
    
    # Insert from temporary table to final table in one query
    cur.execute("""
        INSERT INTO gene (hgnc_symbol, ensembl_gene_id, description, chromosome_name, strand, transcript_count, gene_start, gene_end, unique_exon_count, gene_length)
        SELECT hgnc_symbol, ensembl_gene_id, description, chromosome_name, strand, transcript_count, gene_start, gene_end, unique_exon_count, gene_length
        FROM temp_genes
    """)
    
    inserted_count = cur.rowcount
    
    # Drop temporary table
    cur.execute("DROP TABLE temp_genes")
    
    print(f"Inserted {inserted_count} genes")
    return inserted_count

def insert_rna_consensus(cur, tsv_path: Path = CONS_TSV, replace: bool = False) -> int:
    """ETL: Load RNA consensus TSV into protein_atlas_expression.
    """

    df = load_table(tsv_path, sep='\t')


    # Create temporary load table
    schema = """
        ensembl_gene_id TEXT NOT NULL,
        tissue TEXT NOT NULL,
        nTPM REAL
    """
    
    rows = [
        (row.Gene, row.Tissue, None if pd.isna(row.nTPM) else float(row.nTPM))
        for row in df.itertuples(index=False)
    ]
    
    insert_sql = "INSERT INTO temp_rna_consensus (ensembl_gene_id, tissue, nTPM) VALUES (?, ?, ?)"
    create_temp_table_and_insert(cur, "temp_rna_consensus", schema, rows, insert_sql)

    
    insert_sql = (
        "INSERT INTO protein_atlas_expression (gene_id, tissue, nTPM) "
        "SELECT g.id, t.tissue, t.nTPM "
        "FROM temp_rna_consensus t JOIN gene g ON g.ensembl_gene_id = t.ensembl_gene_id "
        "LEFT JOIN protein_atlas_expression e ON e.gene_id = g.id AND e.tissue = t.tissue "
    )

    cur.execute(insert_sql)
    inserted = cur.rowcount

    cleanup_temp_table(cur, "temp_rna_consensus")
    print(f"Inserted consensus expression rows: {inserted}")
    return inserted

def insert_cdose_map(cur, csv_path: Path = Path("gene_info/CDoseMap.csv")) -> int:
    """ETL: Load CDoseMap CNV data CSV into CDoseMap table.
    """
    df = load_table(csv_path)


    # Rename columns
    rename_map = {
        'CdoseMap_CNV_Type': 'cnv_type',
        'Cytoband': 'cytoband',
        'Size': 'size',
        'CdoseMap_Discovery Sig.': 'discovery_sig',
        'CdoseMap_Known GD': 'known_gd',
        'gnomAD Constrained Genes': 'gnomad_constrained_genes',
        'DB': 'db_source',
        'Article': 'article'
    }
    df = df.rename(columns=rename_map)

    # Temp table
    schema = """
        ensembl_gene_id TEXT,    
        cnv_type TEXT,
        cytoband TEXT,
        size INTEGER,
        discovery_sig TEXT,
        known_gd TEXT,
        gnomad_constrained_genes TEXT,
        db_source TEXT,
        article TEXT
    """

    # Bulk insert raw rows as-is (no per-field None coercion)
    selected_cols = [
        'ensembl_gene_id','cnv_type','cytoband','size',
        'discovery_sig','known_gd','gnomad_constrained_genes','db_source','article'
    ]

    # Ensure all expected columns exist
    check_columns(df, selected_cols)

    rows = list(df[selected_cols].itertuples(index=False, name=None))

    insert_sql = """
        INSERT INTO temp_cdose_map (
            ensembl_gene_id, cnv_type, cytoband, size,
            discovery_sig, known_gd, gnomad_constrained_genes, db_source, article
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    
    create_temp_table_and_insert(cur, "temp_cdose_map", schema, rows, insert_sql)


    insert_sql = """
        INSERT INTO CDoseMap (
            gene_id, cnv_type, cytoband,
            size, discovery_sig, known_gd, gnomad_constrained_genes, db_source, article
        )
        SELECT 
            g.id AS gene_id ,
            t.cnv_type,
            t.cytoband,
            t.size,
            t.discovery_sig,
            t.known_gd,
            t.gnomad_constrained_genes,
            t.db_source,
            t.article
        FROM temp_cdose_map t
        RIGHT JOIN gene g ON g.ensembl_gene_id = t.ensembl_gene_id
        
    """
    cur.execute(insert_sql)
    inserted = cur.rowcount

    cleanup_temp_table(cur, "temp_cdose_map")
    print(f"Inserted CDoseMap rows: {inserted}")
    return inserted



def insert_gtex_gene_stats_enriched(cur, csv_path: Path = Path("gene_info/gtex_gene_stats_enriched.csv")) -> int:
    """ETL: Load GTEx gene statistics enriched data CSV into gtex_gene_stats table.
    """
    df = load_table(csv_path)

    # Source-specific cleaning
    # Remove rows where there is a _PAR_ suffix https://asia.ensembl.org/info/genome/genebuild/human_PARS.html
    df = df[~df['Name'].str.contains('_PAR_')]
    # Clean Ensembl IDs - remove version suffix (e.g., ENSG00000000003.15 -> ENSG00000000003)
    df['ensembl_gene_id'] = df['Name'].str.split('.').str[0]

    # Create temporary table
    schema = """
        ensembl_gene_id TEXT NOT NULL,
        mean REAL,
        median REAL,
        std REAL,
        tau REAL,
        n_zero INTEGER,
        n_below_low INTEGER,
        n_above_low INTEGER,
        n_above_high INTEGER
    """

    # Prepare data for insertion
    selected_cols = ['ensembl_gene_id', 'mean', 'median', 'std', 'tau', 'n_zero', 'n_below_low', 'n_above_low', 'n_above_high']
    # Ensure all expected columns exist
    check_columns(df, selected_cols)
    rows = list(df[selected_cols].itertuples(index=False, name=None))    
    insert_sql = """
        INSERT INTO temp_gtex_stats (
            ensembl_gene_id, mean, median, std, tau,
            n_zero, n_below_low, n_above_low, n_above_high
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    
    create_temp_table_and_insert(cur, "temp_gtex_stats", schema, rows, insert_sql)

    # Insert from temp table to final table, joining with gene table
    insert_sql = """
        INSERT INTO gtex_gene_stats (
            gene_id, mean, median, std, tau,
            n_zero, n_below_low, n_above_low, n_above_high
        )
        SELECT 
            g.id AS gene_id,
            t.mean,
            t.median,
            t.std,
            t.tau,
            t.n_zero,
            t.n_below_low,
            t.n_above_low,
            t.n_above_high
        FROM temp_gtex_stats t
        JOIN gene g ON g.ensembl_gene_id = t.ensembl_gene_id
    """
    
    cur.execute(insert_sql)
    inserted = cur.rowcount

    cleanup_temp_table(cur, "temp_gtex_stats")
    print(f"Inserted GTEx gene stats rows: {inserted}")
    return inserted

def insert_gtex_tpm(cur, csv_path: Path = Path("gene_info/gtex_tissue_tpm.csv")) -> int:
    """ETL: Load GTEx TPM wide table (tissues as columns) into gtex_tissue_tpm.

    Expected header start: Description, Name, <Tissue1>, <Tissue2>, ...
    We treat 'Name' as Ensembl gene ID (strip version) and melt the rest.
    """
    df = load_table(csv_path)


    ensembl_col = df.columns[1]
    tissue_cols = df.columns[2:]
    df = df.rename(columns={ensembl_col: 'ensembl_gene_id'})
    df['ensembl_gene_id'] = df['ensembl_gene_id'].astype(str).str.split('.').str[0]

    long_df = df.melt(id_vars=['ensembl_gene_id'], value_vars=tissue_cols,
                      var_name='tissue', value_name='TPM')
    long_df = long_df[~long_df['TPM'].isna()]

    schema = """
        ensembl_gene_id TEXT NOT NULL,
        tissue TEXT NOT NULL,
        TPM REAL
    """
    rows = [
        (r.ensembl_gene_id, r.tissue, float(r.TPM))
        for r in long_df.itertuples(index=False)
    ]
    insert_sql = "INSERT INTO temp_gtex_tpm (ensembl_gene_id, tissue, TPM) VALUES (?, ?, ?)"
    create_temp_table_and_insert(cur, "temp_gtex_tpm", schema, rows, insert_sql)

    final_sql = """
        INSERT INTO gtex_tissue_tpm (gene_id, tissue, TPM)
        SELECT g.id, t.tissue, t.TPM
        FROM temp_gtex_tpm t
        JOIN gene g ON g.ensembl_gene_id = t.ensembl_gene_id
    """
    cur.execute(final_sql)
    inserted = cur.rowcount
    cleanup_temp_table(cur, "temp_gtex_tpm")
    print(f"Inserted GTEx tissue TPM rows (wide->long): {inserted}")
    return inserted

def insert_drugs_dgi(cur, csv_path: Path = Path("gene_info/drugs_dgi_db.csv")) -> int:
    """ETL: Load drugs DGI aggregated metrics CSV into drugs_dgi table.
    """
    df = load_table(csv_path)

    selected_cols = [
        'gene_name', 'n_drugs', 'n_interaction_types', 'mean_interaction_score',
        'max_interaction_score', 'sum_interaction_score', 'n_approved_drugs',
        'has_immunotherapy', 'has_anti_neoplastic', 'n_types'
    ]
    check_columns(df, selected_cols)

    schema = """
        gene_name TEXT NOT NULL,
        n_drugs INTEGER,
        n_interaction_types INTEGER,
        mean_interaction_score REAL,
        max_interaction_score REAL,
        sum_interaction_score REAL,
        n_approved_drugs INTEGER,
        has_immunotherapy INTEGER,
        has_anti_neoplastic INTEGER,
        n_types INTEGER
    """

    rows = list(df[selected_cols].itertuples(index=False, name=None))
    insert_sql = """
        INSERT INTO temp_drugs_dgi (
            gene_name, n_drugs, n_interaction_types, mean_interaction_score,
            max_interaction_score, sum_interaction_score, n_approved_drugs,
            has_immunotherapy, has_anti_neoplastic, n_types
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    create_temp_table_and_insert(cur, "temp_drugs_dgi", schema, rows, insert_sql)

    insert_final_sql = """
        INSERT INTO drugs_dgi (
            gene_id, n_drugs, n_interaction_types, mean_interaction_score,
            max_interaction_score, sum_interaction_score, n_approved_drugs,
            has_immunotherapy, has_anti_neoplastic, n_types
        )
        SELECT g.id, t.n_drugs, t.n_interaction_types, t.mean_interaction_score,
               t.max_interaction_score, t.sum_interaction_score, t.n_approved_drugs,
               t.has_immunotherapy, t.has_anti_neoplastic, t.n_types
        FROM temp_drugs_dgi t
        JOIN gene g ON g.hgnc_symbol = t.gene_name
    """
    cur.execute(insert_final_sql)
    inserted = cur.rowcount
    cleanup_temp_table(cur, "temp_drugs_dgi")
    print(f"Inserted drugs_dgi rows: {inserted}")
    return inserted


def insert_drugs_central(cur, csv_path: Path = Path("gene_info/drug_central.csv")) -> int:
    """ETL: Load DrugCentral per-gene summary metrics into drugs_central table.
    """
    df = load_table(csv_path)

    # Filter out multi-gene composite entries
    df = df[~df['GENE'].str.contains('\|')]

    # Rename to snake_case matching schema
    rename_map = {
        'GENE': 'gene_name',
        'mean_ACT_VALUE': 'mean_act_value',
        'median_ACT_VALUE': 'median_act_value',
        'n_ACT_TYPE': 'n_act_type',
        'n_TARGET_CLASS': 'n_target_class',
        'best_DRUG_NAME': 'best_drug_name',
        'best_ACT_VALUE': 'best_act_value',
        'best_ACT_TYPE': 'best_act_type',
        'best_TARGET_CLASS': 'best_target_class'
    }
    df = df.rename(columns=rename_map)

    selected_cols = [
        'gene_name', 'n_drugs', 'mean_act_value', 'median_act_value', 'n_act_type',
        'n_target_class', 'best_drug_name', 'best_act_value', 'best_act_type', 'best_target_class'
    ]
    check_columns(df, selected_cols)

    schema = """
        gene_name TEXT NOT NULL,
        n_drugs INTEGER,
        mean_act_value REAL,
        median_act_value REAL,
        n_act_type INTEGER,
        n_target_class INTEGER,
        best_drug_name TEXT,
        best_act_value REAL,
        best_act_type TEXT,
        best_target_class TEXT
    """
    rows = list(df[selected_cols].itertuples(index=False, name=None))
    insert_sql = """
        INSERT INTO temp_drugs_central (
            gene_name, n_drugs, mean_act_value, median_act_value, n_act_type,
            n_target_class, best_drug_name, best_act_value, best_act_type, best_target_class
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    create_temp_table_and_insert(cur, "temp_drugs_central", schema, rows, insert_sql)

    final_insert_sql = """
        INSERT INTO drugs_central (
            gene_id, n_drugs, mean_act_value, median_act_value, n_act_type,
            n_target_class, best_drug_name, best_act_value, best_act_type, best_target_class
        )
        SELECT g.id, t.n_drugs, t.mean_act_value, t.median_act_value, t.n_act_type,
               t.n_target_class, t.best_drug_name, t.best_act_value, t.best_act_type, t.best_target_class
        FROM temp_drugs_central t
        JOIN gene g ON g.hgnc_symbol = t.gene_name
    """
    cur.execute(final_insert_sql)
    inserted = cur.rowcount
    cleanup_temp_table(cur, "temp_drugs_central")
    print(f"Inserted drugs_central rows: {inserted}")
    return inserted

def insert_gnomad(cur, csv_path: Path = Path("gene_info/gnoMAD_full_corrected.csv")) -> int:
    """ETL: Load gnomAD constraint metrics (v2.1.1 & v4.1) into gnomad table.

    Source columns (original):
        mis.oe_v2.1.1, mis.mu_v2.1.1, mis_pphen.oe_v2.1.1, lof.mu_v2.1.1, lof.pLI_v2.1.1,
        lof.oe_v2.1.1, mis.z_score_v2.1.1, lof.z_score_v2.1.1, gene_id,
        lof.oe_v4.1, lof.pLI_v4.1, lof.z_raw_v4.1, mis.oe_v4.1, mis.z_raw_v4.1, mis_pphen.oe_v4.1, ensembl_gene_id

    We keep both gene_id and ensembl_gene_id columns, but join using ensembl_gene_id
    after stripping version (if present). Duplicate ensembl IDs are collapsed by
    choosing the first occurrence (dataset appears unique already).
    """
    df = load_table(csv_path)
    
    #check that values in df['ensembl_gene_id'] and df['gene_id']  are the same, then drop df['gene_id']
    if not df['ensembl_gene_id'].equals(df['gene_id']):
        raise ValueError("Ensembl IDs and Gene IDs do not match")
    non_unique = df[df.duplicated(subset=['ensembl_gene_id'], keep=False)] #84 rows
    #removing non-unique rows, TODO this needs an investigation
    df = df[~df.duplicated(subset=['ensembl_gene_id'], keep=False)]
    df = df.drop(columns=['gene_id'])

    rename_map = {
        'mis.oe_v2.1.1': 'mis_oe_v2_1_1',
        'mis.mu_v2.1.1': 'mis_mu_v2_1_1',
        'mis_pphen.oe_v2.1.1': 'mis_pphen_oe_v2_1_1',
        'lof.mu_v2.1.1': 'lof_mu_v2_1_1',
        'lof.pLI_v2.1.1': 'lof_pLI_v2_1_1',
        'lof.oe_v2.1.1': 'lof_oe_v2_1_1',
        'mis.z_score_v2.1.1': 'mis_z_score_v2_1_1',
        'lof.z_score_v2.1.1': 'lof_z_score_v2_1_1',
        'lof.oe_v4.1': 'lof_oe_v4_1',
        'lof.pLI_v4.1': 'lof_pLI_v4_1',
        'lof.z_raw_v4.1': 'lof_z_raw_v4_1',
        'mis.oe_v4.1': 'mis_oe_v4_1',
        'mis.z_raw_v4.1': 'mis_z_raw_v4_1',
        'mis_pphen.oe_v4.1': 'mis_pphen_oe_v4_1'
    }
    df = df.rename(columns=rename_map)

    selected_cols = [
        'ensembl_gene_id', 'mis_oe_v2_1_1', 'mis_mu_v2_1_1', 'mis_pphen_oe_v2_1_1',
        'lof_mu_v2_1_1', 'lof_pLI_v2_1_1', 'lof_oe_v2_1_1', 'mis_z_score_v2_1_1', 'lof_z_score_v2_1_1',
        'lof_oe_v4_1', 'lof_pLI_v4_1', 'lof_z_raw_v4_1', 'mis_oe_v4_1', 'mis_z_raw_v4_1', 'mis_pphen_oe_v4_1'
    ]
    check_columns(df, selected_cols)

    schema = """
        ensembl_gene_id TEXT NOT NULL,
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
        mis_pphen_oe_v4_1 REAL
    """

    rows = [tuple(row) for row in df[selected_cols].itertuples(index=False, name=None)]
    insert_sql = """
        INSERT INTO temp_gnomad (
            ensembl_gene_id, mis_oe_v2_1_1, mis_mu_v2_1_1, mis_pphen_oe_v2_1_1,
            lof_mu_v2_1_1, lof_pLI_v2_1_1, lof_oe_v2_1_1, mis_z_score_v2_1_1, lof_z_score_v2_1_1,
            lof_oe_v4_1, lof_pLI_v4_1, lof_z_raw_v4_1, mis_oe_v4_1, mis_z_raw_v4_1, mis_pphen_oe_v4_1
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    create_temp_table_and_insert(cur, "temp_gnomad", schema, rows, insert_sql)

    final_insert = """
        INSERT INTO gnomad (
            gene_id, mis_oe_v2_1_1, mis_mu_v2_1_1, mis_pphen_oe_v2_1_1,
            lof_mu_v2_1_1, lof_pLI_v2_1_1, lof_oe_v2_1_1, mis_z_score_v2_1_1, lof_z_score_v2_1_1,
            lof_oe_v4_1, lof_pLI_v4_1, lof_z_raw_v4_1, mis_oe_v4_1, mis_z_raw_v4_1, mis_pphen_oe_v4_1
        )
        SELECT g.id, t.mis_oe_v2_1_1, t.mis_mu_v2_1_1, t.mis_pphen_oe_v2_1_1,
               t.lof_mu_v2_1_1, t.lof_pLI_v2_1_1, t.lof_oe_v2_1_1, t.mis_z_score_v2_1_1, t.lof_z_score_v2_1_1,
               t.lof_oe_v4_1, t.lof_pLI_v4_1, t.lof_z_raw_v4_1, t.mis_oe_v4_1, t.mis_z_raw_v4_1, t.mis_pphen_oe_v4_1
        FROM temp_gnomad t
        JOIN gene g ON g.ensembl_gene_id = t.ensembl_gene_id
    """
    cur.execute(final_insert)
    inserted = cur.rowcount
    cleanup_temp_table(cur, "temp_gnomad")
    print(f"Inserted gnomad rows: {inserted}")
    return inserted


def main():
    """Main migration function."""
    print("Starting migration from CSV/TSV to SQLite...")
    print(f"Database will be created at: {DB_PATH}")

    conn, cur = create_database()
    try:
        migrate_genes(cur)
        insert_rna_consensus(cur)
        insert_cdose_map(cur)
        insert_gtex_gene_stats_enriched(cur)
        insert_gtex_tpm(cur)
        insert_drugs_dgi(cur)
        insert_drugs_central(cur)
        insert_gnomad(cur)
        conn.commit()
        print("Migration completed successfully")
    except Exception as exc:
        print(f"Migration failed: {exc}")
        conn.rollback()
        raise
    finally:
        conn.close()

if __name__ == "__main__":
    main()
