import pandas as pd
import sys
import sqlite3
from pathlib import Path
DB_PATH = "genes.db"
SCHEMA_PATH = Path("schema/schema.sql")


def check_columns(df, selected_cols):
    for col in selected_cols:
        if col not in df.columns:
            print(f"Column '{col}' not found in input data")
            sys.exit()

def create_temp_table_and_insert(cur, table_name: str, schema: str, data: list, insert_sql: str) -> None:
    """Generic temporary table creation and data insertion.
    
    Args:
        cur: SQLite cursor
        table_name: Name of the temporary table to create
        schema: SQL schema definition for the table columns
        data: List of tuples containing the data to insert
        insert_sql: SQL INSERT statement with placeholders
    """
    cur.execute(f"DROP TABLE IF EXISTS {table_name}")
    cur.execute(f"CREATE TEMPORARY TABLE {table_name} ({schema})")
    cur.executemany(insert_sql, data)

def cleanup_temp_table(cur, table_name: str) -> None:
    """Drop temporary table.
    
    Args:
        cur: SQLite cursor
        table_name: Name of the temporary table to drop
    """
    cur.execute(f"DROP TABLE {table_name}")

def load_table(file_path: Path, sep: str = ',') -> pd.DataFrame:
    """Load CSV/TSV with standardized error handling.
    
    Args:
        file_path: Path to the CSV/TSV file
        sep: Separator character (default: ',')
        
    Returns:
        DataFrame with loaded data, or empty DataFrame if file doesn't exist or fails to load
    """
    if not file_path.exists():
        print(f"WARNING: File not found at {file_path}")
        sys.exit()
    try:
        return pd.read_csv(file_path, sep=sep)
    except Exception as e:
        print(f"ERROR: Failed to read {file_path}: {e}")
        sys.exit()

def create_database():
    """Create the SQLite database with schema."""
    print("Creating database and schema...")
    
    if not SCHEMA_PATH.exists():
        print(f"ERROR: Schema file not found at {SCHEMA_PATH}")
        sys.exit(1)
    
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    
    # Load and execute schema
    with open(SCHEMA_PATH, "r", encoding="utf-8") as f:
        schema_sql = f.read()
    
    cur.executescript(schema_sql)
    conn.commit()
    print("✓ Database schema created")
    
    return conn, cur