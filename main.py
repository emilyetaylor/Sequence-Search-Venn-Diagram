# Author: Emily Taylor, for the Data Science Team in the Goetsch Lab at Michigan Technological University
# Email: eetaylor@mtu.edu
# Purpose: This main.py file serves as the driver script for the sequence search venn diagram pipeline, automating
# searches to blastx, tblastn, and hmmer, combining the resulting hits, and creating a venn diagram of the resulting
# hits based on their search method of origin.

import blast_query
import hmmer_query
import pandas as pd
import json
import map_accessions

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # example consensus sequence:
    consensus_species = input("Enter model species: ")
    protein = input("Enter consensus sequence protein name: ")
    accession = input("Enter consensus sequence: ")
    # Start blastp search

    record = blast_query.fetch_protein_sequence(accession)
    sequence = blast_query.convert_seqrecord_to_fasta(record)
    print(sequence)
    # print(consensus_species + " " + protein + " record: " + record)
    # result_handle_blastp = blast_query.run_blastp(record)
    # blastp_df =blast_query.parse_blast_results(result_handle_blastp)
    #
    # print(blastp_df.head())
    #
    # # Start blastx search
    # result_handle_blastx = blast_query.run_blastx(record)
    # blastx_df = blast_query.parse_blast_results(result_handle_blastx)
    #
    # print(blastx_df.head())

    # Start hmmer search
    job_id = hmmer_query.submit_hmmer_search(sequence)
    results = hmmer_query.wait_for_completion(job_id)
    print(json.dumps(results, indent=2))

    hmmer_df = hmmer_query.parse_results(results)

    print(hmmer_df.head())

    # Get proper ids from uniprot ID mapping service

    job_id = map_accessions.submit_id_mapping("UniProtKB_AC-ID","RefSeq_Protein", hmmer_df)
    results_json = map_accessions.get_results(job_id)
    df = map_accessions.parse_mapped_results(results_json)
    print(df.head())

    # Start comparison of accessions and venn diagram
    accession_df = ...      # should be titled using name of consensus species, protein, and original accession
    vennDiagram_img = ...   # should be titled using name of consensus species, protein, and original accession

