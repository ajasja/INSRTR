import pandas as pd
import pickle


def encode_categories(df, replace=False):
    """
    Encodes categories using cat.codes from scikit.

    Parameters
    ----------
    df: input dataframe
    replace: if True, the original categorical features will be dropped from the dataset

    Returns
    -------
    modified dataframe
    """
    cols_to_encode = [
        "resi_type",
        "resi_dssp",
        "prev_resi_type",
        "prev_resi_dssp",
        "next_resi_type",
        "next_resi_dssp",
        "loop_seq",
    ]
    encoded_cols = [col + "_encoded" for col in cols_to_encode]
    # Apply astype and cat.codes to each column
    encoded = df[cols_to_encode].apply(lambda x: x.astype("category").cat.codes)
    # Use assign to create new columns in dataframe
    df = df.assign(**dict(zip(encoded_cols, encoded.T.values)))
    # Drop original columns if specified
    if replace:
        df.drop(cols_to_encode, axis=1, inplace=True)
    return df


def load_model(filename):
    """
    Loads a scikit-learn model from disk using the pickle module.

    Parameters
    ----------
    filename: where the pickle file is saved

    Returns
    -------
    loaded model
    """
    with open(filename, "rb") as file:
        model = pickle.load(file)
    return model


def predict_positions(df, model_path="models/gbt_classifier_v1.pkl", n_top=3):
    """
    Loads the model and applies it to the input dataframe.
    Takes into account only the positive predictions.
    Only the highest probability prediction for each loop_index0 is considered.

    Parameters
    ----------
    df: input dataframe
    model_path: path to trained model
    n_top: the number of to positions to return

    Returns
    -------
    df_predictions: dataframe with top 3 predictions
    df_all
    """
    # Load the model
    model = load_model(model_path)
    # Preprocess the data
    df_preprocess = df.drop(columns=["struct_name"])
    x = encode_categories(df_preprocess, replace=True).values
    # Apply model to get labels and predicted probabilities
    prediction_label = model.predict(x)
    prediction_probability = model.predict_proba(x)
    # Add probabilities for N and Y to dataframe
    df["probability_N"] = prediction_probability[:, 0]
    df["probability_Y"] = prediction_probability[:, 1]
    # Make a dataframe with only positive predictions
    df_positive = df[prediction_label == "Y"].rename(columns={"probability_Y": "prediction_probability"})
    # Create the output dataframes
    df_predictions = df_positive.loc[
        df_positive.groupby("loop_index0")["prediction_probability"].idxmax().sample(frac=1, random_state=2),
        ["resi_index0", "resi_dssp", "prediction_probability"],
    ].nlargest(n=n_top, columns=["prediction_probability"])
    return df_predictions, df

