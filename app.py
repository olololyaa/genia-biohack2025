import streamlit as st
import pandas as pd
import sqlite3
from typing import Optional
import numpy as np 
from plots import create_violin_plot

st.set_page_config(page_title="GIA", layout="centered", page_icon="🧬",)
st.title("Gene Insights & Analysis")

st.markdown(
    """
    <style>
    [data-testid="stSidebar"] > div:first-child {
        display: flex;
        flex-direction: column;
        align-items: center;
        justify-content: center;
    }
    </style>
    """,
    unsafe_allow_html=True
)

@st.cache_resource(show_spinner=False)
def get_conn() -> Optional[sqlite3.Connection]:
    try:
        # Allow usage across Streamlit threads
        conn = sqlite3.connect("genes.db", check_same_thread=False)
        conn.row_factory = sqlite3.Row
        return conn
    except Exception as e:
        st.error(f"Could not connect to SQLite database: {e}")
        return None

@st.cache_data(show_spinner=False)
def load_gene_table() -> Optional[pd.DataFrame]:
    conn = get_conn()
    if not conn:
        return None
    try:
        df = pd.read_sql_query(
            "SELECT hgnc_symbol, ensembl_gene_id, description, chromosome_name FROM gene", conn
        )
        return df
    except Exception as e:
        st.error(f"Failed to load gene data from database: {e}")
        return None

@st.cache_data(show_spinner=False)
def load_consensus_for_ensembl(ensembl_ids: list[str]) -> Optional[pd.DataFrame]:
    if not ensembl_ids:
        return pd.DataFrame(columns=["Gene", "Tissue", "nTPM"])
    conn = get_conn()
    if not conn:
        return None
    try:
        placeholders = ",".join(["?"] * len(ensembl_ids))
        query = f"""
            SELECT g.ensembl_gene_id AS Gene, e.tissue AS Tissue, e.nTPM
            FROM protein_atlas_expression e
            JOIN gene g ON g.id = e.gene_id
            WHERE g.ensembl_gene_id IN ({placeholders})
        """
        df = pd.read_sql_query(query, conn, params=ensembl_ids)
        return df
    except Exception as e:
        st.error(f"Failed to load consensus data: {e}")
        return None

# --- New: CDoseMap loader ---
@st.cache_data(show_spinner=False)
def load_cdosemap_for_ensembl(ensembl_ids: list[str]) -> Optional[pd.DataFrame]:
    if not ensembl_ids:
        return pd.DataFrame(columns=["cnv_type","cytoband","size","discovery_sig","known_gd","gnomad_constrained_genes","db_source","article"])    
    conn = get_conn()
    if not conn:
        return None
    try:
        placeholders = ",".join(["?"] * len(ensembl_ids))
        q = f"""
            SELECT 
                   c.cnv_type,
                   c.cytoband,
                   c.size,
                   c.discovery_sig,
                   c.known_gd,
                   c.gnomad_constrained_genes,
                   c.db_source,
                   c.article
            FROM CDoseMap c
            JOIN gene g ON g.id = c.gene_id
            WHERE g.ensembl_gene_id IN ({placeholders})
        """
        return pd.read_sql_query(q, conn, params=ensembl_ids).dropna().drop_duplicates()
    except Exception as e:
        st.error(f"Failed to load CDoseMap data: {e}")
        return None



# --- New: Gtex loader ---
@st.cache_data(show_spinner=False)
def load_gtex_data(ensembl_ids: list[str] = None) -> Optional[pd.DataFrame]:
    if not ensembl_ids:
        query = f"""
            SELECT g.ensembl_gene_id, mean, median, std, tau, n_zero, n_below_low, n_above_low, n_above_high
            FROM gtex_gene_stats e
            JOIN gene g ON g.id = e.gene_id
        """
    else:
        placeholders = ",".join(["?"] * len(ensembl_ids))
        query += "WHERE g.ensembl_gene_id IN ({placeholders})"
    conn = get_conn()
    if not conn:
        return None
    try:
        df = pd.read_sql_query(query, conn, params=ensembl_ids)
        return df
    except Exception as e:
        st.error(f"Failed to load GTEx data: {e}")
        return None

# Replace top-level UI with tabbed layout
gene_tab, about_tab = st.tabs(["Gene Search", "About"])

with gene_tab:
    df = load_gene_table()

    if df is not None and not df.empty:
        with st.sidebar:
            st.subheader("Gene Search")
            user_input = st.text_input(
                "Enter genes (HGNC symbols or Ensembl IDs); separate multiple entries with commas",
                placeholder="e.g., BRCA1 or ENSG00000141510",
                key="manual_input",
            )
            selected_genes = user_input.replace(" ", "").split(",") if user_input else []
            # Validate terms: keep only those present in either hgnc_symbol or ensembl_gene_id.
            # Collect those not found for user feedback.
            if selected_genes:
                df_symbols = set(df['hgnc_symbol'].dropna().astype(str)) if 'hgnc_symbol' in df.columns else set()
                df_ensembl = set(df['ensembl_gene_id'].dropna().astype(str)) if 'ensembl_gene_id' in df.columns else set()
                valid_set = df_symbols | df_ensembl
                # Preserve original order while filtering
                filtered = [t for t in selected_genes if t and t in valid_set]
                not_found = [t for t in selected_genes if t and t not in valid_set]
                selected_genes = filtered
                if not selected_genes:
                    st.warning("No valid gene identifiers found.")
                elif not_found:
                    st.warning(f"Not found / invalid: {', '.join(not_found)}")

        if not selected_genes:
            st.markdown("Please select at least one gene identifier to view details.")
        else:
            st.subheader(f"Selected genes {', '.join(selected_genes)}")
            display_cols = ['hgnc_symbol','ensembl_gene_id','chromosome_name','description']
            
            sel_rows = df[df['hgnc_symbol'].isin(selected_genes) | df['ensembl_gene_id'].isin(selected_genes)]
            existing = [c for c in display_cols if c in sel_rows.columns]
            for _, row in sel_rows[existing].iterrows():
                title = row.get('hgnc_symbol') or row.get('ensembl_gene_id')
                st.markdown(f"### {title}")
                for c in existing:
                    st.markdown(f"**{c}**: {row[c] if pd.notna(row[c]) else ''}")
                st.markdown("---")

            sel_ensembl_ids = sel_rows['ensembl_gene_id'].dropna().astype(str).str.strip().unique().tolist() if 'ensembl_gene_id' in sel_rows.columns else []

            with st.expander("Consensus RNA data"):
                consensus_df = load_consensus_for_ensembl(sel_ensembl_ids)
                if consensus_df is None or consensus_df.empty:
                    st.info("No consensus expression data found for the selected gene(s).")
                else:
                    st.markdown("[Consensus RNA data from the Human Protein Atlas](https://www.proteinatlas.org/humanproteome/tissue/data#consensus_tissues_rna)")
                    st.markdown("The consensus normalized expression ('nTPM') value is calculated as the maximum nTPM value for each gene in the two data sources.")
                    st.dataframe(consensus_df, use_container_width=True, hide_index=True)

            with st.expander("CDoseMap data"):
                cd_df = load_cdosemap_for_ensembl(sel_ensembl_ids)
                if cd_df is None or cd_df.empty:
                    st.markdown("No CDoseMap data found for the selected gene(s).")
                    st.image("tg_image_3445653099.jpeg", width=200)
                else:
                    show_cols = [c for c in [
                        'cnv_type','cytoband','size','discovery_sig','known_gd',
                        'gnomad_constrained_genes','db_source','article'
                    ] if c in cd_df.columns]
                    cd_clean = cd_df.dropna(how='all', subset=show_cols)
                    if cd_clean.empty:
                        st.info("No CNV annotations (only empty rows).")
                    else:
                        st.markdown("Data from [the paper](https://www.cell.com/cell/fulltext/S0092-8674(22)00788-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422007887%3Fshowall%3Dtrue#mmc1)")
                        st.markdown("Data were provided using gene symbols only; do not extrapolate to specific Ensembl IDs.")
                        st.dataframe(cd_clean[show_cols], use_container_width=True, hide_index=True)

            with st.expander("Enriched Gene Statistics"): #TODO naming!
                gene_symbols_for_plot = sel_rows['ensembl_gene_id']

                if len(gene_symbols_for_plot) == 0:
                    st.info("No Ensembl IDs available for the selected entries.")
                else:
                    # read gene_info\gene_stats_enriched.csv
                    #enriched_gene_stats = pd.read_csv("gene_info/gtex_gene_stats_enriched.csv")
                    enriched_gene_stats = load_gtex_data()
                    #import pdb; pdb.set_trace()
                    #ensemble version ids are provided within the Name column
                    # split values in the name column by "." and leave only first element
                    #enriched_gene_stats['Name'] = enriched_gene_stats['Name'].str.split(".").str[0]
                    st.caption("Sample synthetic data (not real expression).")
                    st.caption("Two plots will be shown: one for the selected gene(s) and another for the selected feature group.")
                    features = ["Feature 1", "Feature 2"] # TODO: define features from other datasets, for example is this gene is drug target or not
                    selected_feature = st.selectbox("Select a feature group", options=features)
                    st.caption("Note: Data for the selected genes will be excluded from the second plot, even if they share the same feature.")
                    stats = ["mean", "median", "std", "tau", "n_zero", "n_below_low", "n_above_low", "n_above_high"]
                    for stat in stats:
                        create_violin_plot(gene_symbols_for_plot, plot_df=enriched_gene_stats, stat=stat, features=features)

            with st.expander("Selected Genes on Сhromosomes"):
                # Upload gene_info\db_paralogues_ortholog_mouse.csv
                dtypes = {i: 'str' for i in range(13)}
                db_paralogues_path = Path("gene_info/db_paralogues_ortholog_mouse.csv")
                db_paralogues_mouse = pd.read_csv(db_paralogues_path, dtype=dtypes, usecols=[0, 7])
                db_paralogues_mouse.drop_duplicates('Gene stable ID', inplace=True)
                # get list of ensembl_gene_ids based on selected_genes list:
                selected_genes_ens = sel_rows['ensembl_gene_id']
                plot_count_genes_on_chromosomes(selected_genes_ens, db_paralogues_mouse)


with about_tab:
    try:
        with open("APP_README.md", "r", encoding="utf-8") as f:
            content = f.read()
        st.markdown(content)
    except FileNotFoundError:
        st.warning("APP_README.md not found in the project root.")
    except Exception as e:
        st.error(f"Could not read APP_README.md: {e}")


