# Author: Emily Taylor, for the Data Science Team in the Goetsch Lab at Michigan Technological University
# Email: eetaylor@mtu.edu
# Purpose: This hmmer_query.py file contains all the functions pertaining to performing a hmmer search utilizing REST
# APIs, querying the ebi server. The end product is a pandas dataframe containing all information pertinent to the
# search that was extracted from the results json.

import requests
import time
import pandas as pd
import json

# submit_hmmer_search(fasta_sequence) - This function will submit a fasta file sequence and a string denoting the
# database to be used.
# INPUT: fasta_sequence - the fasta file sequence record to be used as the query
# RETURN: resp.json()["id"] -  the job ID, usually a string or int
def submit_hmmer_search(fasta_sequence):
    url = "https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/phmmer"
    headers = {
        "User-Agent": "PythonScript/1.0",
        "Accept" : "application/json",
        "Content-Type" : "application/json"
    }
    payload = {
        "input": fasta_sequence,
        "database" : "refprot"  # option, default is 'uniprotrefprot' and can be changed
    }

    resp = requests.post(url, json=payload, headers=headers)
    resp.raise_for_status()
    return resp.json()["id"]

# wait_for_completion(job_id) - this function will continually check the status of the job being performed
# on the ebi server and will display a spinner while the job is running. This works by continually submitting requests
# to the status url every three seconds. An error will be raised if an unknown status is returned.
# INPUT job_id: a string containing the unique job id from the submission (performed in submit_hmmer_search)
# RETURN result_resp.json() - the json object containing the results from the initial query. Recieved from the result_url
def wait_for_completion(job_id, max_retries = 10):
    status_url = f"https://www.ebi.ac.uk/Tools/hmmer/api/v1/result/{job_id}"
    result_url = f"https://www.ebi.ac.uk/Tools/hmmer/api/v1/result/{job_id}"
    headers = {"Accept" : "application/json"}

    retries = 0
    while True:
        resp = requests.get(status_url, headers=headers)
        resp.raise_for_status()
        data = resp.json()
        status = data.get("status", "")

        if status in ("STARTED", "RUNNING", "STARTED"):
            print("Still running...")
            time.sleep(3)
        elif status == "RETRY":
            if retries < max_retries:
                retries+=1
                print(f"Status is RETRY - attempt {retries}/{max_retries}")
                print("Raw status response:", json.dumps(data, indent=2))
                time.sleep(5)
            else:
                raise RuntimeError("Exceeded max retries for RETRY status.")
        elif status == "SUCCESS":
            print("Job completed. Fetching results...")
            result_resp = requests.get(result_url, headers = headers)
            result_resp.raise_for_status()
            return result_resp.json()
        else:
            raise RuntimeError(f"Job failed or returned status: {status}")


# parse_results(results_json) - this function will extract relevant fields of data from the json and parse the data into
# a pandas dataframe. Specific data is specified and can be changed according to pipeline needs.
# INPUT results_json - the json file object returned from the hmmer query. Contains all data
# RETURN df - pandas dataframe containing specified data from json
def parse_results(results_json):
    # Go into "result" → "hits"
    result = results_json.get("result", {})
    hits = result.get("hits", [])

    print(f"Number of hits: {len(hits)}")
    if hits:
        print("First hit example:")
        print(hits[0])  # show first full hit
        print("Metadata of first hit:")
        print(hits[0].get("metadata", {}))

    # Extract relevant fields from the "metadata" of each hit
    parsed_data = []
    for hit in hits:
        metadata = hit.get("metadata", {})
        row = {     # types of data that are extracted can be changed/modified according to pipeline needs
            "accession": metadata.get("accession"),
            "identifier": metadata.get("identifier"),
            "uniprot accession": metadata.get("uniprot_accession"),
            "phylum": metadata.get("phylum"),
            "species": metadata.get("species"),
        }
        parsed_data.append(row)
    # Convert to DataFrame
    df = pd.DataFrame(parsed_data)
    return df