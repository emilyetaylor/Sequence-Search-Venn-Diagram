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

# Test sequence: C. elegans lin-9 NP_001023016

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # example consensus sequence:
    consensus_species = input("Enter model species: ")
    protein = input("Enter consensus sequence protein name: ")
    accession = input("Enter consensus sequence: ")
    # Start blastp search

    record = blast_query.fetch_protein_sequence(accession)
    sequence = blast_query.convert_seqrecord_to_fasta(record)
    # result_handle_blastp = blast_query.run_blastp(record)
    # blastp_df =blast_query.parse_blast_results(result_handle_blastp)
    #
    # print(blastp_df.head())
    #
    # # Start blastx search
    # result_handle_blastx = blast_query.run_blastx(record)
    # blastx_df = blast_query.parse_blast_results(result_handle_blastx)




    # Start hmmer search
    job_id = hmmer_query.submit_hmmer_search(sequence)
    results = hmmer_query.wait_for_completion(job_id)

    hmmer_df = hmmer_query.parse_results(results)

    print(hmmer_df.head())

    # Get proper ids from uniprot ID mapping service

# Start hmmer mapping
    job_id = hmmer_query.submit_hmmer_search(sequence)
    results = hmmer_query.wait_for_completion(job_id)
    hmmer_df = hmmer_query.parse_results(results)

    print(hmmer_df.head())

    # Get proper ids from uniprot ID mapping service

    #Hmmer mapping
    job_id = map_accessions.submit_id_mapping("UniProtKB_AC-ID","RefSeq_Protein", hmmer_df)
    results_json = map_accessions.get_results(job_id)
    df, failed_ids = map_accessions.parse_mapped_results(results_json)
    print(f"Mapped Ids: {len(df)}, Failed Ids: {len(failed_ids)}")
    print(f"Failed Ids: {failed_ids}")

    #hmmer_results = df + failed_ids

    # Blastp

    blastp_results, blastp_failed_ids = map_accessions.blast_mapping(blastp_df)

    # Blastx

    blastx_results, blastx_failed_ids = map_accessions.blast_mapping(blastx_df)


    # Start comparison of accessions and venn diagram

    vennDiagram_img = ...   # should be titled using name of consensus species, protein, and original accession