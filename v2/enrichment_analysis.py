# This script takes a list of gene symbols as input and performs enrichment analysis
# against Gene Ontology (GO), KEGG, and Reactome pathway databases.

import gseapy as gp
import pandas as pd
import logging

# Suppress verbose logging from gseapy to keep the output clean
logging.basicConfig(level=logging.WARNING)

def perform_enrichment_analysis(gene_list):
    """
    Performs enrichment analysis on a given gene list using gseapy's enrichr function.

    Args:
        gene_list (list): A list of gene symbols (e.g., ['STAT3', 'IL6', 'TNF']).

    Returns:
        None. Prints the results and saves them to CSV files.
    """
    if not gene_list:
        print("The gene list is empty. Please provide a list of genes.")
        return

    print("Starting enrichment analysis for {} genes...".format(len(gene_list)))

    # Define the databases you want to query.
    # You can find more libraries using: gp.get_library_name()
    gene_sets = [
        'GO_Biological_Process_2023',
        'GO_Cellular_Component_2023',
        'GO_Molecular_Function_2023',
        'KEGG_2021_Human',
        'Reactome_2022'
    ]

    try:
        # Run the enrichment analysis using Enrichr
        # organism='Human' is the default, but it's good to be explicit.
        # Other options include 'Mouse', 'Yeast', 'Fly', 'Fish', 'Worm'.
        enr = gp.enrichr(gene_list=gene_list,
                         gene_sets=gene_sets,
                         organism='Human',
                         outdir='enrichment_results', # Directory to save the results
                         cutoff=0.05, # P-value cutoff
                         no_plot=True # We will generate a custom plot later
                        )

        print("\nEnrichment analysis complete. Results are saved in the 'enrichment_results' directory.")

        # The result object `enr` has an attribute `results` which is a pandas DataFrame.
        # Let's inspect the results for one of the databases.
        for db in gene_sets:
            print("\n--- Top 10 results for {} ---".format(db))
            
            # Access results for each database by its name
            db_results_df = enr.results[enr.results['Gene_set'] == db]
            
            if not db_results_df.empty:
                # Print the top 10 significant terms
                print(db_results_df.head(10))
                
                # You can save the full results for each database to a separate CSV
                # The files are already saved by gseapy in the 'enrichment_results' folder,
                # but this shows how you could do it manually if needed.
                # db_results_df.to_csv('{}_results.csv'.format(db), index=False)

            else:
                print("No significant terms found for {} with the current settings.".format(db))


        # --- Visualization ---
        # Create a dot plot for the top 10 results from a specific database.
        # Dot plots are great for visualizing term, p-value, and gene count simultaneously.
        print("\nGenerating dot plot for 'GO_Biological_Process_2023'...")
        
        # Check if the key exists and has results
        go_bp_results = enr.results[enr.results['Gene_set'] == 'GO_Biological_Process_2023']
        if not go_bp_results.empty:
            # You can also filter the results before plotting
            # For example, plot only terms with adjusted p-value < 0.05
            significant_go_bp = go_bp_results[go_bp_results['Adjusted P-value'] < 0.05]

            if not significant_go_bp.empty:
                # Let's create a custom title for the plot
                plot_title = "Top Enriched GO Biological Processes"

                # The `dotplot` function can take the results DataFrame directly
                ax = gp.dotplot(significant_go_bp.head(15), # Plot top 15 significant terms
                              title=plot_title,
                              cutoff=0.05,
                              ofname='enrichment_dotplot.png') # Save the plot
                print("Dot plot saved as 'enrichment_dotplot.png'")
            else:
                print("No significant GO Biological Process terms to plot.")

    except Exception as e:
        print("An error occurred during the analysis: {}".format(e))
        print("This may be due to an issue with the gene list or connection to the Enrichr servers.")


if __name__ == '__main__':
    # --- INPUT YOUR GENE LIST HERE ---
    # This is a sample gene list related to the IL-6/JAK/STAT3 signaling pathway.
    # Replace it with your own list of gene symbols.
    my_gene_list = [
        'STAT3', 'IL6', 'TNF', 'IL1B', 'CXCL8', 'VEGFA', 'EGF', 'EGFR', 'JAK1', 'JAK2',
        'SOCS3', 'MAPK1', 'MAPK3', 'HIF1A', 'NFKB1', 'RELA', 'JUN', 'FOS', 'MYC', 'BCL2L1'
    ]

    perform_enrichment_analysis(my_gene_list)
