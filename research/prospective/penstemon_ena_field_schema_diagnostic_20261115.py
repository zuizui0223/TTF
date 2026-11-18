"""Diagnose public ENA returnFields response envelope ONLY (no biological records).

The first two attempts to parse ENA 'result=sample' returnFields as JSON
or ordinary TSV did not match the actual response. This source-only
diagnosis records MIME, body length, and a bounded first-line prefix,
not sequence, sample or phenotype records.
"""
import json
import urllib.error
import urllib.request

URL="https://www.ebi.ac.uk/ena/portal/api/returnFields?result=sample"
LIMIT=2_000_000


def main():
    req=urllib.request.Request(URL,headers={
        "User-Agent":"TTF-Habroanthus-ENA-field-schema-diagnosis/1.0",
        "Accept":"*/*",
    })
    try:
        with urllib.request.urlopen(req,timeout=35) as resp:
            data=resp.read(LIMIT+1)
            code=resp.status
            ctype=resp.headers.get("Content-Type")
    except urllib.error.HTTPError as exc:
        print(json.dumps({"status":"SCHEMA_ENDPOINT_HTTP_ERROR",
              "http_status":exc.code,"no_sample_data_retrieved":True}))
        return
    if len(data)>LIMIT:raise ValueError("metadata dictionary exceeded safety cap")
    value=data.decode("utf-8-sig",errors="replace")
    # The endpoint is a published dictionary, not sample data.
    line1=value.splitlines()[0][:250] if value.splitlines() else ""
    line2=value.splitlines()[1][:250] if len(value.splitlines())>1 else ""
    print(json.dumps({
        "status":"ENA_PUBLIC_FIELD_SCHEMA_RESPONSE_DIAGNOSIS_ONLY",
        "url":URL,"http_status":code,"mime":ctype,
        "response_bytes":len(data),"line_count":len(value.splitlines()),
        "first_schema_line_prefix":line1,
        "second_schema_line_prefix":line2,
        "first_char_code_point":ord(value[0]) if value else None,
        "no_biological_data_retrieved":True
    },indent=2,sort_keys=True))


if __name__=="__main__":main()
