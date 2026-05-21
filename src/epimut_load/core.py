from pathlib import Path
import pandas as pd
import numpy as np
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



def filter_meta(meta, include_filter="", exclude_filter=""):

    meta_filt = meta.copy()

    print("Metadata dimensions: ", meta_filt.shape)

    if include_filter != "":
        print("Inclusion filter:", include_filter)
        meta_filt = meta_filt.query(include_filter, engine="python")

    if exclude_filter != "":
        print("Exclusion filter:", exclude_filter, engine="python")
        meta_filt = meta_filt.query(f"not ({exclude_filter})")

    if include_filter != "" or exclude_filter != "":
        print("Filtered metadata dimensions: ", meta_filt.shape)

    return meta_filt



def select_ref_samples(meta, config):

    sample_id_col = config["input"]["sample_id_col"]
    mode = config["reference"]["mode"]

    if mode == "internal":

        include_filter = config["reference"]["internal"].get("include_filter", "")
        exclude_filter = config["reference"]["internal"].get("exclude_filter", "")

        meta_ref = filter_meta(
            meta,
            include_filter=include_filter,
            exclude_filter=exclude_filter
        )

    else:
        raise ValueError("External reference is not implemented yet.")

    ref_samples = meta_ref[sample_id_col].tolist()

    if len(ref_samples) == 0:
        raise ValueError("No reference samples selected.")

    return ref_samples



def select_target_samples(meta, config):

    sample_id_col = config["input"]["sample_id_col"]

    include_filter = config["target_samples"].get("include_filter", "")
    exclude_filter = config["target_samples"].get("exclude_filter", "")

    meta_target = filter_meta(
        meta,
        include_filter=include_filter,
        exclude_filter=exclude_filter
    )

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



def build_annots(outliers):

    d_so = {}

    for col in outliers.columns:

        d_so[col] = {}

        epimut = outliers[outliers[col] != 0].index.tolist()
        hyper = outliers[outliers[col] == 1].index.tolist()
        hypo = outliers[outliers[col] == -1].index.tolist()

        d_so[col]["epimut"] = epimut
        d_so[col]["hyper"] = hyper
        d_so[col]["hypo"] = hypo

    return d_so



def compute_eml(meta, outliers, sem_strength, config):

    sample_id_col = config["input"]["sample_id_col"]

    eml = outliers.abs().sum()
    sum_sem = outliers.sum()

    meta_out = meta[meta[sample_id_col].isin(outliers.columns)].copy()

    meta_out["eml"] = meta_out[sample_id_col].map(eml)
    meta_out["log1p_eml"] = meta_out["eml"].apply(lambda x: np.log1p(x))
    meta_out["sum_sem"] = meta_out[sample_id_col].map(sum_sem)

    if sem_strength is not None:

        sum_pos_strength = sem_strength[sem_strength > 0].sum(skipna=True)
        sum_neg_strength = sem_strength[sem_strength < 0].sum(skipna=True)

        meta_out["sum_pos_strength"] = meta_out[sample_id_col].map(sum_pos_strength)
        meta_out["sum_neg_strength"] = meta_out[sample_id_col].map(sum_neg_strength)

    return meta_out



def save_outputs(outliers, sem_strength, d_so, meta_out, config):

    output_path = Path(config["output"]["path"])
    name_exp = config["output"]["name_exp"]

    output_path.mkdir(parents=True, exist_ok=True)

    outliers.to_csv(output_path / f"outliers_{name_exp}.tsv", sep="\t")

    if sem_strength is not None:
        sem_strength.to_csv(output_path / f"sem_strength_{name_exp}.tsv", sep="\t")

    pd.DataFrame(d_so).transpose().to_csv(
        output_path / f"sems_{name_exp}.tsv",
        sep="\t"
    )

    meta_out.to_csv(
        output_path / f"sample_file_metrics_{name_exp}.tsv",
        sep="\t",
        index=False
    )



### MAIN
def run_pipeline(config_path):

    # read config file
    config = load_config(config_path)

    # load data
    print("\nLoading input files...")
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


    # reference
    print("\n\n********************")
    print("REFERENCE: Selecting samples...")
    print("********************")
    ref_samples = select_ref_samples(meta, config)
    # print(f"Reference samples: {len(ref_samples)}")

    # target
    print("\n\n********************")
    print("TARGET: Selecting samples...")
    print("********************")
    target_samples = select_target_samples(meta, config)
    # print(f"Target samples: {len(target_samples)}")

    
    # compute ref quantiles
    print("\n\nComputing reference quantiles...")
    quantiles = get_quantiles(data, ref_samples)
    print(quantiles)

    # detect SEMs
    print("\n\nDetecting SEMs...")
    iqr_multiplier = config["sem_calling"]["iqr_multiplier"]

    data_target = data[target_samples].copy()

    outliers = data_target.apply(
        detect_outliers,
        axis=1,
        quantiles=quantiles,
        iqr_multiplier=iqr_multiplier
    )

    print("\nOutliers matrix:")
    print(outliers.shape)
    print(outliers.iloc[:5, :5])


    # compute SEMs strength
    sem_strength = None

    if config["sem_calling"].get("compute_strength", True):

        print("\n\nComputing SEM strength...")

        sem_strength = data_target.apply(
            compute_sem_strength,
            axis=1,
            quantiles=quantiles,
            iqr_multiplier=iqr_multiplier
        )

        print("\nSEM strength matrix:")
        print(sem_strength.shape)
        print(sem_strength.iloc[:5, :5])



    print("\n\nBuilding SEM annotations...")
    d_so = build_annots(outliers)

    print("\n\nComputing EML...")
    meta_out = compute_eml(meta, outliers, sem_strength, config)

    print("\n\nSaving outputs...")
    save_outputs(outliers, sem_strength, d_so, meta_out, config)

    print("\nDone.")