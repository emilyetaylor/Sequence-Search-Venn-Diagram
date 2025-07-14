# Author: Emily Taylor, for the Data Science Team in the Goetsch Lab at Michigan Technological University
# Email: eetaylor@mtu.edu
# Purpose: This map_accessions.py file contains all the functions pertaining to performing an id mapping to a different
# database utilizing REST APIs, querying the uniprot server. The end product is a pandas dataframe containing all
# mapping information pertinent to making sure that results from hmmer, blastx, and tblastn are comparable.

import requests
import pandas as pd
import time
# submit_id_mapping(from_db, to_db, df) - this function will submit an list of accessions, to uniprot, to get each accession's
# unique id for another database. It will return the number of ids submitted for mapping, and will return a job id string
# INPUT from_db - the database where the provided ids are coming from
# INPUT to_db - the database where the provided ids need to be mapped to (given new ids)
# INPUT df - the dataframe containing the original accessions to be mapped
# RETURN response.json()["jobID"] - job id of the mapping done by uniprot
def submit_id_mapping(from_db, to_db, df):
    accession_list = df["accession"].dropna().unique().tolist()
    url = "https://rest.uniprot.org/idmapping/run"
    headers = {
        "Content-Type" : "application/x-www-form-urlencoded"
    }
    payload = {
        "from" : from_db,
        "to" : to_db,
        "ids" : ",".join(accession_list)
    }

    print("=== SUBMITTING ID MAPPING REQUEST ===")
    print("From:", from_db)
    print("To:", to_db)
    print("Total IDs submitted:", len(accession_list))
    response = requests.post(url, data = payload, headers=headers)
    response.raise_for_status()
    return response.json()["jobId"]

# get_entry(job_id) - this function checks the status of the job on the server side and updates via the console.
# INPUT(job_id) - the id of the specific job, taken from the json returned by the initial post to the uniprot server to
# initialize the query
# RETURN void - this function simply stalls the execution of further code until the job is finished
def get_entry(job_id):
    while True:
        status = check_status(job_id)
        if status["jobStatus"] == "FINISHED":
            break
        print("Waiting for mapping to finish")
        time.sleep(3)

# check_status(job_id) - this function directly checks the status of the job from the job id. This function is used in
# the get_entry function and is continually called until the job is finished.
# INPUT job_id - the job id (string) used to access the status url
# RETURN response.json() - the json object that contains the information about the status of the job.
# Returns to get_entry function as that is sole usage
# TODO: Implement check_status recursively utilizing code from get_entry. Get rid of get_entry
def check_status(job_id):
    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    response = requests.get(status_url)
    response.raise_for_status()
    return response.json()

# get_results(job_id) - this function returns a json containing the results of the given job id
# INPUT job_id - the specific job id provided by the submit_id_mapping function
# RETURN response.json() - the json object containing the results of the query
def get_results(job_id):
    result_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
    response = requests.get(result_url)
    response.raise_for_status()
    return response.json()

# parse_mapped_results(json_data) - creates a pandas dataframe containing all the mapping information of the mapped ids
# (database from & associated id, database to & associated id)
# INPUT json_data - the mapping data in the json, provided by the get_results function
# RETURN pd.DataFrame(parsed) - a pandas datafram containing the parsed information from the json (see description above)
def parse_mapped_results(json_data):
    results = json_data.get("results", [])
    failed_ids = json_data.get("failedIds", [])
    parsed = []

    for entry in results:
        from_id = entry["from"]
        to_id = entry["to"]
        parsed.append({"uniprot": from_id, "mapped_id": to_id})


    return pd.DataFrame(parsed), failed_ids

def blast_mapping(df):
    blast_mapped = []     # refseq
    blast_mapping = []
    for entry in df:
        if entry['hit_id'].startswith('ref'):
            blast_mapping.append(entry)
        elif entry['hit_id'].startswith('gb', 'emb', 'dbj'):
            blast_mapping.append(entry)

    job_id = submit_id_mapping("EMBL-GenBank-DDBJ", "RefSeq_Protein", blast_mapping)

    results_json = get_results(job_id)
    blast_mapped_df, blast_failed_ids = parse_mapped_results(results_json)

    blast_mapped_df+=blast_mapped

    return blast_mapped_df, blast_failed_ids


