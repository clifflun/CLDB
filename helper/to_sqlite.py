"""
Load TSV/CSV files into SQLite with configurable engine, inputs, metadata, and table names.
Edit CONFIG below to swap datasets and run multiple reload jobs.
"""
import glob
import yaml
import pandas as pd
import sqlite3 as sq
from sqlalchemy import create_engine
from dataclasses import dataclass, field
from typing import Optional, List

pd.set_option('display.expand_frame_repr', False)


# ---------------------------------------------------------------------------
# Config: swap engine, input files, metadata file, table names here
# ---------------------------------------------------------------------------

@dataclass
class TableLoadJob:
    """One TSV → SQLite table load (supports chunked reading)."""
    input_glob: str
    table_name: str
    chunksize: int = 100_000
    index_columns: Optional[list[str]] = None  # columns to create indexes on


@dataclass
class MetaLoadJob:
    """One metadata TSV → SQLite table (full read)."""
    input_file: str
    table_name: str
    columns: list[str]  # e.g. ['pt_id', 'sex', 'phenotype']


@dataclass
class Config:
    engine_url: str 
    table_jobs: list[TableLoadJob] = field(default_factory=list)
    meta_jobs: list[MetaLoadJob] = field(default_factory=list)


def load_config_from_yaml(filepath: str) -> List[Config]:
    """Hook to read YAML and transform it into a list of Config objects."""
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)
    
    configs = []
    for db_entry in data.get('databases', []):
        # Build TableLoadJobs
        table_jobs = [
            TableLoadJob(**job) for job in db_entry.get('table_jobs', [])
        ]
        
        # Build MetaLoadJobs
        meta_jobs = [
            MetaLoadJob(**job) for job in db_entry.get('meta_jobs', [])
        ]
        
        # Create the main Config object
        configs.append(Config(
            engine_url=f"sqlite:///{db_entry['db_name']}",
            table_jobs=table_jobs,
            meta_jobs=meta_jobs
        ))
    return configs


# ---------------------------------------------------------------------------
# Core functions (all take engine/params; no hardcoded paths)
# ---------------------------------------------------------------------------

def get_engine(engine_url: str):
    return create_engine(engine_url)


def write_meta_to_DB(engine_url: str, input_file: str, table_name: str, columns: list[str]):
    df = pd.read_csv(input_file, sep='\t')
    df = df[columns]
    print(df)
    engine = get_engine(engine_url)
    df.to_sql(table_name, con=engine, if_exists='replace', index=False)


def write_to_DB(
    engine_url: str,
    input_glob: str,
    table_name: str,
    chunksize: int = 100_000,
):
    engine = get_engine(engine_url)
    fnames = glob.glob(input_glob)
    print(fnames)
    for fn in fnames:
        print(fn)
        reader = pd.read_csv(fn, sep='\t', chunksize=chunksize, low_memory=False)
        for idx, chunk in enumerate(reader):
            print(f'processing chunk {idx}')
            chunk.to_sql(table_name, con=engine, if_exists='append', index=False)


def create_index(engine_url: str, table_name: str, index_columns: list[str]):
    # SQLite URL is "sqlite:///path" -> path is after 3 slashes
    db_path = engine_url.replace('sqlite:///', '')
    conn = sq.connect(db_path)
    cursor = conn.cursor()
    for col in index_columns:
        print(col)
        cursor.execute(f'CREATE INDEX IF NOT EXISTS {col}_index ON {table_name} ({col})')
    conn.commit()
    conn.close()


def delete_table(engine_url: str, table_name: str):
    db_path = engine_url.replace('sqlite:///', '')
    print(f'deleting table {table_name}')
    conn = sq.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(f'DROP TABLE IF EXISTS {table_name}')
    conn.commit()
    conn.close()
    print('deleted table')


def run_config(config: Config, drop_tables_first: bool = True):
    """Run all table and metadata jobs from config."""
    for job in config.table_jobs:
        if drop_tables_first:
            delete_table(config.engine_url, job.table_name)
        write_to_DB(
            config.engine_url,
            job.input_glob,
            job.table_name,
            job.chunksize,
        )
        if job.index_columns:
            create_index(config.engine_url, job.table_name, job.index_columns)

    for job in config.meta_jobs:
        write_meta_to_DB(
            config.engine_url,
            job.input_file,
            job.table_name,
            job.columns,
        )


def main():
    # Now you just point to your external config file
    config_file = 'Z:/Members/clun/CLDB/helper/sqlite_config.yaml'
    
    try:
        configs = load_config_from_yaml(config_file)
        
        for config in configs:
            print(f"--- Starting Load: {config.engine_url} ---")
            run_config(config, drop_tables_first=True)
            
    except FileNotFoundError:
        print(f"Error: {config_file} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    main()
