import json
import pandas as pd
import os


# Define the list of tools that will become columns in the DataFrame
TOOLS_LIST = ["miRanda", "miRmap", "RNAhybrid", "PITA", "TargetScan"]

def process_json_to_dataframe(json_file_path):
    """
    Reads a JSON file, extracts data from the "prediction" key, and converts it
    into a Pandas DataFrame.

    The "prediction" key in the JSON should be a dictionary where keys are
    tool names (e.g., "miRanda", "PITA") and values are lists of gene names
    predicted by that tool.

    Args:
        json_file_path (str): The path to the input JSON file.

    Returns:
        pd.DataFrame | None: A DataFrame with gene names and tool predictions,
                              or None if an error occurs.
    """
    try:
        with open(json_file_path, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print("Error: The file {} was not found.".format(json_file_path))
        return None
    except json.JSONDecodeError:
        print("Error: Could not decode JSON from the file {}.".format(json_file_path))
        return None
    except Exception as e:
        print("An unexpected error occurred while reading the file: {}".format(e))
        return None

    predictions_data = data.get("prediction")
    if predictions_data is None:
        print("Error: 'prediction' key not found in {}.".format(json_file_path))
        return None
    if not isinstance(predictions_data, dict):
        print("Error: 'prediction' key in {} is not a dictionary (expected tool names as keys).".format(json_file_path))
        return None

    # This dictionary will store gene_name -> set_of_tools_that_predicted_it
    gene_to_tools_map = {}

    for tool_name_from_json, predicted_genes_list in predictions_data.items():
        # We only consider tools that are defined in our TOOLS_LIST
        if tool_name_from_json not in TOOLS_LIST:
            # Optionally, print a warning if a tool in JSON is not in TOOLS_LIST
            # print("Warning: Tool '{}' from JSON is not in the recognized TOOLS_LIST and will be ignored.".format(tool_name_from_json))
            continue

        if not isinstance(predicted_genes_list, list):
            print("Warning: Skipping tool '{}' as its prediction data (gene list) is not a list.".format(tool_name_from_json))
            continue
        
        for gene_name in predicted_genes_list:
            if not isinstance(gene_name, str):
                print("Warning: Skipping non-string gene identifier '{}' found under tool '{}'.".format(gene_name, tool_name_from_json))
                continue
            gene_to_tools_map.setdefault(gene_name, set()).add(tool_name_from_json)

    if not gene_to_tools_map:
        print("No valid gene prediction data found to process after parsing tools and genes.")
        return pd.DataFrame(columns=["gene_name"] + TOOLS_LIST)

    processed_data_list = []
    all_gene_names = sorted(list(gene_to_tools_map.keys())) # Sort for consistent output order

    for gene_name in all_gene_names:
        row = {"ENST_ID": gene_name}
        predicting_tools_for_this_gene = gene_to_tools_map.get(gene_name, set())
        for tool_in_df_column in TOOLS_LIST:
            row[tool_in_df_column] = 1 if tool_in_df_column in predicting_tools_for_this_gene else 0
        processed_data_list.append(row)
    
    df_columns = ["ENST_ID"] + TOOLS_LIST
    df = pd.DataFrame(processed_data_list, columns=df_columns)
    
    return df

def process_dataframe(df, coding_map):
    df = coding_map.merge(df, on='ENST_ID', how='left')
    df = df.dropna()
    df["SUM"] = df['miRanda'] + df['PITA'] + df['miRmap'] + df['RNAhybrid'] + df['TargetScan']
    return df

def main():
    path = ''
    for file in os.listdir(path):
        file_path = os.path.join(path, file)
        df = process_json_to_dataframe(file_path)
        coding_map = pd.read_csv('coding_map.tsv', sep='\t')
        #non_coding_map = pd.read_csv('non_coding_map.tsv', sep='\t')
        df = process_dataframe(df, coding_map)
        output_path = os.path.join(path, file.replace('.json', '.tsv'))
        if df is not None:
            if not df.empty:
                try:
                    df.to_csv(output_path, index=False, sep='\t')
                    print("\nDataFrame successfully saved to {}".format(file.replace('.json', '.tsv')))
                except Exception as e:
                    print("\nError saving DataFrame to CSV: {}".format(e))

if __name__ == "__main__":
    main() 
