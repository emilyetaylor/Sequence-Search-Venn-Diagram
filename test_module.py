# Author: Emily Taylor, for the Data Science Team in the Goetsch Lab at Michigan Technological University
# Email: eetaylor@mtu.edu
# Purpose: This test_module.py file contains all the functions pertaining to performing all automated blastx, tblastn,
# and hmmer queries, as well as extracting information from the query's respective jsons and xml results files, and
# parsing the resulting data into pandas dataframes. The resulting data will then be mapped so the ids of all hit
# sequences are comparable, which will serve to combine all hit sequences from all search methods and help to create
# a venn diagram showcasing the differing results of each search method.

import pytest
import pandas as pd
from unittest.mock import patch, Mock, MagicMock
from hmmer_query import submit_hmmer_search, wait_for_completion, parse_results
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from blast_query import (
    fetch_protein_sequence,
    convert_seqrecord_to_fasta,
    run_blastp,
    get_mrna_from_protein,
    parse_blast_results
)
from map_accessions import (
    submit_id_mapping,
    get_entry,
    check_status,
    get_results,
    parse_mapped_results
)

# -------------------------------
# Mock data for hmmer_query tests
# -------------------------------

FASTA_SEQUENCE = ">example\nMTEITAAMVKELRESTGAGMMDCKNALSETQHEVLHGLTD"

MOCK_JOB_ID = "job123"

MOCK_STATUS_RESPONSE_STARTED = {
    "status": "STARTED"
}

MOCK_STATUS_RESPONSE_SUCCESS = {
    "status": "SUCCESS"
}

MOCK_RESULT_RESPONSE = {
    "result": {
        "hits": [
            {
                "metadata": {
                    "accession": "P12345",
                    "identifier": "SAMPLE_PROT",
                    "uniprot_accession": "UP12345",
                    "phylum": "Proteobacteria",
                    "species": "Escherichia coli"
                }
            }
        ]
    }
}

# -------------------------------
# hmmer_query.py tests
# -------------------------------

# This function tests the submit_hmmer_search function in hmmer_query.py.
# It passes when the correct job_id is returned given a mock API post request (quality parameter)
@patch("hmmer_query.requests.post")
def test_submit_hmmer_search(mock_post):
    mock_response = Mock()
    mock_response.json.return_value = {"id": MOCK_JOB_ID}
    mock_response.raise_for_status = Mock()
    mock_post.return_value = mock_response

    job_id = submit_hmmer_search(FASTA_SEQUENCE)
    assert job_id == MOCK_JOB_ID
    mock_post.assert_called_once()

# TODO: create test cases for submit_hmmer_search with bad parameter(s) and edge cases. Analyze errors that are raised.

# This function tests the wait_for_completion function in hmmer_query.py.
# It passes when the first call to the mock api returns 'STARTED', and the second call returns 'SUCCESS', and when the
# result response is returned as an instance, as expected in a specific format.
@patch("hmmer_query.requests.get")
def test_wait_for_completion(mock_get):
    # First call returns STARTED, second returns SUCCESS
    mock_get.side_effect = [
        Mock(status_code=200, json=lambda: MOCK_STATUS_RESPONSE_STARTED, raise_for_status=Mock()),
        Mock(status_code=200, json=lambda: MOCK_STATUS_RESPONSE_SUCCESS, raise_for_status=Mock()),
        Mock(status_code=200, json=lambda: MOCK_RESULT_RESPONSE, raise_for_status=Mock())
    ]

    result = wait_for_completion(MOCK_JOB_ID)
    assert "result" in result
    assert isinstance(result["result"]["hits"], list)

# TODO: create test cases for wait_for_completion with bad parameter(s) and edge cases. Analyze errors that are raised.

# This function tests the parse_results function in hmmer_query.py.
# It passes when the format of the data from the mock api request matches the specified format.
@patch("builtins.print")  # to suppress print output during test
def test_parse_results(mock_print):
    df = parse_results(MOCK_RESULT_RESPONSE)

    assert isinstance(df, pd.DataFrame)
    assert df.shape[0] == 1
    assert df.iloc[0]["accession"] == "P12345"
    assert df.iloc[0]["species"] == "Escherichia coli"

# TODO: create test cases for parse_results with bad inputs and edge cases. Analyze errors that are raised.

# -------------------------------
# Mock data for blast_query tests
# -------------------------------
# Mock protein sequence record to be used for testing blast_query.py
mock_protein = SeqRecord(Seq("MKTIIALSYIFCLVFADYKDDDDK"), id="NP_000240.1", description="Mock protein")

# -----------------------------
# blast_query.py tests
# -----------------------------

# This function tests the fetch_protein_sequence function in blast_query.py.
# It passes when the correct protein sequence ncbi accession number is returned.
@patch("blast_query.Entrez.efetch")
@patch("blast_query.SeqIO.read")
def test_fetch_protein_sequence(mock_read, mock_efetch):
    mock_read.return_value = mock_protein
    mock_efetch.return_value = MagicMock()

    result = fetch_protein_sequence("NP_000240.1")

    assert result.id == "NP_000240.1"
    mock_efetch.assert_called_once()

# TODO: create test cases for fetch_protein_sequence with bad inputs and edge cases. Analyze errors that are raised.

# This function tests the convert_seqrecord_to_fasta function in blast_query.py.
# It passes when the a string containing the correct protein sequence is returned.
def test_convert_seqrecord_to_fasta():
    fasta_str = convert_seqrecord_to_fasta(mock_protein)
    assert fasta_str.startswith(">NP_000240.1")
    assert "MKTIIALSYIFCLVFADYKDDDDK" in fasta_str


# TODO: create test cases for convert_seqrecord_to_fasta with bad inputs and edge cases. Analyze errors that are raised.

# This function tests the run_blastp function in blast_query.py.
# It passes when the result returned from the mock API is not None - e.g. the function simply needs to return something.
# TODO: edit this test to increase specificity in terms of result. What should the function actually be returning?
@patch("blast_query.NCBIWWW.qblast")
@patch("blast_query.yaspin")  # Suppress spinner for test
def test_run_blastp(mock_spinner, mock_qblast):
    mock_spinner.return_value.__enter__.return_value = MagicMock(ok=MagicMock())
    mock_qblast.return_value = MagicMock()

    result = run_blastp(mock_protein)
    assert result is not None
    mock_qblast.assert_called_once()

# TODO: create test cases for run_blastp with bad inputs and edge cases. Analyze errors that are raised.

# This function tests the get_mrna_from_protein function in blast_query.py.
# It passes when the result id returned from the mock api matches the expected result id.
@patch("blast_query.SeqIO.read")
@patch("blast_query.Entrez.efetch")
@patch("blast_query.Entrez.elink")
@patch("blast_query.Entrez.read")
def test_get_mrna_from_protein(mock_read, mock_elink, mock_efetch, mock_seqread):
    # Set up mock link response
    mock_read.return_value = [{"LinkSetDb": [{"Link": [{"Id": "123456"}]}]}]
    mock_elink.return_value = MagicMock()
    mock_seqread.return_value = SeqRecord(Seq("ATGC"), id="XM_123456.1", description="Mock mRNA")
    mock_efetch.return_value = MagicMock()

    result = get_mrna_from_protein(mock_protein)
    assert result.id == "XM_123456.1"

# TODO: create more test cases for get_mrna_from_protein with bad inputs and edge cases. Analyze errors that are raised.

# -----------------------------
# parse_blast_results
# -----------------------------
@patch("blast_query.NCBIXML.read")
def test_parse_blast_results(mock_read):
    # Mock HSP and alignment
    hsp = MagicMock()
    hsp.expect = 0.001
    hsp.identities = 50
    hsp.align_length = 100
    hsp.bits = 200
    hsp.query_start = 1
    hsp.query_end = 100
    hsp.sbjct_start = 1
    hsp.sbjct_end = 100

    alignment = MagicMock()
    alignment.hit_id = "sp|P12345|"
    alignment.hit_def = "Mock protein"
    alignment.accession = "P12345"
    alignment.hsps = [hsp]

    mock_blast_record = MagicMock()
    mock_blast_record.query = "MockQuery"
    mock_blast_record.alignments = [alignment]

    mock_read.return_value = mock_blast_record

    from io import StringIO
    handle = StringIO("mock")

    df = parse_blast_results(handle)
    assert not df.empty
    assert df.iloc[0]["hit_id"] == "sp|P12345|"


# Sample DataFrame for testing
@pytest.fixture
def sample_df():
    return pd.DataFrame({"accession": ["P12345", "Q67890", "A1B2C3", None]})

# -----------------------------
# submit_id_mapping
# -----------------------------
@patch("map_accessions.requests.post")
def test_submit_id_mapping(mock_post, sample_df):
    mock_response = MagicMock()
    mock_response.json.return_value = {"jobId": "mock_job_123"}
    mock_response.raise_for_status = MagicMock()
    mock_post.return_value = mock_response

    job_id = submit_id_mapping("UniProtKB_AC", "RefSeq_Protein", sample_df)
    assert job_id == "mock_job_123"
    mock_post.assert_called_once()

# -----------------------------
# check_status
# -----------------------------
@patch("map_accessions.requests.get")
def test_check_status(mock_get):
    mock_response = MagicMock()
    mock_response.json.return_value = {"jobStatus": "RUNNING"}
    mock_response.raise_for_status = MagicMock()
    mock_get.return_value = mock_response

    status = check_status("mock_job_123")
    assert status["jobStatus"] == "RUNNING"
    mock_get.assert_called_once()

# -----------------------------
# get_entry
# -----------------------------
@patch("map_accessions.time.sleep", return_value=None)
@patch("map_accessions.check_status")
def test_get_entry(mock_check_status, mock_sleep):
    # First call returns not done, second returns done
    mock_check_status.side_effect = [
        {"jobStatus": "RUNNING"},
        {"jobStatus": "FINISHED"}
    ]

    get_entry("mock_job_123")
    assert mock_check_status.call_count == 2
    mock_sleep.assert_called_once()

# -----------------------------
# get_results
# -----------------------------
@patch("map_accessions.requests.get")
def test_get_results(mock_get):
    mock_response = MagicMock()
    mock_response.json.return_value = {"results": [{"from": "P12345", "to": "XP_123456"}]}
    mock_response.raise_for_status = MagicMock()
    mock_get.return_value = mock_response

    results = get_results("mock_job_123")
    assert "results" in results
    assert results["results"][0]["from"] == "P12345"

# -----------------------------
# parse_mapped_results
# -----------------------------
def test_parse_mapped_results():
    sample_json = {
        "results": [
            {"from": "P12345", "to": "XP_123456"},
            {"from": "Q67890", "to": "XP_789012"},
        ]
    }

    df = parse_mapped_results(sample_json)
    assert isinstance(df, pd.DataFrame)
    assert df.shape[0] == 2
    assert set(df.columns) == {"uniprot", "mapped_id"}
    assert df.iloc[0]["mapped_id"] == "XP_123456"