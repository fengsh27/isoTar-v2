# This script takes a list of gene symbols as input and performs enrichment analysis
# against Gene Ontology (GO), KEGG, and Reactome pathway databases.

import os
import gseapy as gp
import pandas as pd
import logging

# Suppress verbose logging from gseapy to keep the output clean
logging.basicConfig(level=logging.WARNING)

def perform_enrichment_analysis(gene_list, organism='Human', cutoff=0.05, outdir='enrichment_results'):
    """
    Performs enrichment analysis on a given gene list using gseapy's enrichr function.

    Args:
        gene_list (list): A list of gene symbols (e.g., ['STAT3', 'IL6', 'TNF']).
        organism (str):   Organism name passed to Enrichr. Options include 'Human',
                          'Mouse', 'Yeast', 'Fly', 'Fish', 'Worm'. Default: 'Human'.
        cutoff (float):   Adjusted p-value cutoff for filtering results. Default: 0.05.
        outdir (str):     Directory to save result CSV files and dot plot. Default: 'enrichment_results'.

    Returns:
        None. Prints the results and saves them to CSV files.
    """
    if not gene_list:
        print("The gene list is empty. Please provide a list of genes.")
        return

    print("Starting enrichment analysis for {} genes (organism={}, cutoff={})...".format(
        len(gene_list), organism, cutoff))

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
        enr = gp.enrichr(gene_list=gene_list,
                         gene_sets=gene_sets,
                         organism=organism,
                         outdir=outdir,
                         cutoff=cutoff,
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
            significant_go_bp = go_bp_results[go_bp_results['Adjusted P-value'] < cutoff]

            if not significant_go_bp.empty:
                # Let's create a custom title for the plot
                plot_title = "Top Enriched GO Biological Processes"

                # The `dotplot` function can take the results DataFrame directly
                dotplot_path = os.path.join(outdir, 'enrichment_dotplot.png')
                ax = gp.dotplot(significant_go_bp.head(15), # Plot top 15 significant terms
                              title=plot_title,
                              cutoff=cutoff,
                              ofname=dotplot_path)
                print("Dot plot saved as '{}'".format(dotplot_path))
            else:
                print("No significant GO Biological Process terms to plot.")

    except Exception as e:
        print("An error occurred during the analysis: {}".format(e))
        print("This may be due to an issue with the gene list or connection to the Enrichr servers.")


if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description='Perform gene set enrichment analysis using Enrichr.'
    )

    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('-f', '--file', type=str,
                        help='Path to a CSV file containing a "gene_label" column')
    source.add_argument('-g', '--genes', type=str, nargs='+',
                        help='Gene symbols to analyse (space-separated)')

    parser.add_argument('-o', '--organism', type=str, default='Human',
                        help='Organism for Enrichr (e.g. Human, Mouse, Fly, Fish, Worm). Default: Human')
    parser.add_argument('-c', '--cutoff', type=float, default=0.05,
                        help='Adjusted p-value cutoff for filtering results. Default: 0.05')
    parser.add_argument('-d', '--outdir', type=str, default='enrichment_results',
                        help='Output directory for result CSV files and dot plot. Default: enrichment_results')

    args = parser.parse_args()

    if args.file:
        df = pd.read_csv(args.file)
        if 'gene_label' not in df.columns:
            print("Error: CSV file must contain a 'gene_label' column.")
            sys.exit(1)
        gene_list = df['gene_label'].dropna().astype(str).str.strip()
        gene_list = [g for g in gene_list if g]
    else:
        gene_list = args.genes

    perform_enrichment_analysis(gene_list, organism=args.organism, cutoff=args.cutoff, outdir=args.outdir)
