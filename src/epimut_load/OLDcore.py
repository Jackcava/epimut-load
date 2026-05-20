from pathlib import Path

import pandas as pd
import yaml


def load_config(config_path):
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    return config


def read_file(path):
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"File not found: {path.resolve()}")

    if path.suffix == ".csv":
        return pd.read_csv(path)

    if path.suffix in [".tsv", ".txt"]:
        return pd.read_csv(path, sep="\t")

    raise ValueError(f"Unsupported file extension: {path.suffix}")



def load_data(config):
    input_path = Path(config["input"]["path"])

    methylation_path = input_path / config["input"]["methylation_matrix"]
    metadata_path = input_path / config["input"]["metadata"]

    methylation = read_file(methylation_path)
    metadata = read_file(metadata_path)

    sample_id_col = config["input"]["sample_id_col"]
    cpg_id_col = config["data"]["cpg_id_col"]

    if sample_id_col not in metadata.columns:
        metadata = metadata.rename(columns={metadata.columns[0]: sample_id_col})

    if cpg_id_col not in methylation.columns:
        methylation = methylation.rename(columns={methylation.columns[0]: cpg_id_col})

    methylation = methylation.set_index(cpg_id_col)

    metadata[sample_id_col] = metadata[sample_id_col].astype(str)
    methylation.columns = methylation.columns.astype(str)

    return methylation, metadata



def select_ref_samples(meta, config):

    sample_id_col = config["input"]["sample_id_col"]
    mode = config["reference"]["mode"]

    if mode == "metadata_query":
        query = config["reference"]["query"]
        meta_ref = meta.query(query)

    elif mode == "all":
        meta_ref = meta.copy()

    else:
        raise ValueError(f"Unsupported reference mode: {mode}")

    ref_samples = meta_ref[sample_id_col].tolist()

    if len(ref_samples) == 0:
        raise ValueError("No reference samples selected.")

    return ref_samples



def select_target_samples(meta, config):

    sample_id_col = config["input"]["sample_id_col"]
    mode = config["target_samples"]["mode"]

    if mode == "all":
        meta_target = meta.copy()

    elif mode == "all_except_missing":
        missing_col = config["target_samples"]["missing_column"]
        meta_target = meta[meta[missing_col].notna()].copy()

    elif mode == "metadata_query":
        query = config["target_samples"]["query"]
        meta_target = meta.query(query)

    else:
        raise ValueError(f"Unsupported target samples mode: {mode}")

    target_samples = meta_target[sample_id_col].tolist()

    if len(target_samples) == 0:
        raise ValueError("No target samples selected.")

    return target_samples



def get_quantiles(data, ref_samples):

    data_ref = data[ref_samples].copy()

    quantiles = data_ref.quantile([0.25, 0.5, 0.75], axis=1)
    quantiles.loc["iqr"] = quantiles.loc[0.75] - quantiles.loc[0.25]

    quantiles = quantiles.transpose()
    quantiles.columns = ["q1", "q2", "q3", "iqr"]

    return quantiles



def compute_sem_strength(row, quantiles, iqr_multiplier):

    median = quantiles.loc[row.name, "q2"]

    lower_bound = quantiles.loc[row.name, "q1"] - iqr_multiplier * quantiles.loc[row.name, "iqr"]
    upper_bound = quantiles.loc[row.name, "q3"] + iqr_multiplier * quantiles.loc[row.name, "iqr"]

    def get_strength(x):
        if x < lower_bound:
            return x - median
        elif x > upper_bound:
            return x - median
        else:
            return 0

    return row.apply(get_strength)



def detect_outliers(row, quantiles, iqr_multiplier):

    lower_bound = quantiles.loc[row.name, "q1"] - iqr_multiplier * quantiles.loc[row.name, "iqr"]
    upper_bound = quantiles.loc[row.name, "q3"] + iqr_multiplier * quantiles.loc[row.name, "iqr"]

    return row.apply(lambda x: -1 if x < lower_bound else (1 if x > upper_bound else 0))



def run_pipeline(config_path):

    config = load_config(config_path)

    print("Loading input files...")
    data, meta = load_data(config)

    sample_id_col = config["input"]["sample_id_col"]

    print("Checking sample matching...")

    meta_samples = set(meta[sample_id_col])
    data_samples = set(data.columns)

    common_samples = meta_samples.intersection(data_samples)

    print(f"Samples in metadata: {len(meta_samples)}")
    print(f"Samples in methylation matrix: {len(data_samples)}")
    print(f"Common samples: {len(common_samples)}")

    if len(common_samples) == 0:
        raise ValueError("No matching samples between metadata and methylation matrix.")

    print("Selecting reference samples...")
    ref_samples = select_ref_samples(meta, config)

    print("Selecting target samples...")
    target_samples = select_target_samples(meta, config)

    print(f"Reference samples: {len(ref_samples)}")
    print(f"Target samples: {len(target_samples)}")

    print("Computing reference quantiles...")
    quantiles = get_quantiles(data, ref_samples)

    print("Detecting SEMs...")
    iqr_multiplier = config["sem_calling"]["iqr_multiplier"]

    data_target = data[target_samples].copy()

    outliers = data_target.apply(
        detect_outliers,
        axis=1,
        quantiles=quantiles,
        iqr_multiplier=iqr_multiplier
    )

    print("Outliers matrix:")
    print(outliers.shape)
    print(outliers.iloc[:5, :5])

    sem_strength = None

    if config["sem_calling"].get("compute_strength", True):

        print("Computing SEM strength...")

        sem_strength = data_target.apply(
            compute_sem_strength,
            axis=1,
            quantiles=quantiles,
            iqr_multiplier=iqr_multiplier
        )

        print("SEM strength matrix:")
        print(sem_strength.shape)
        print(sem_strength.iloc[:5, :5])

    print("\nDone.")