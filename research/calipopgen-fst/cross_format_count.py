"""Compare original CSV/XLSX structural row counts without interpreting GenDist."""
import argparse,csv,hashlib,io,json,urllib.request
from pathlib import Path
from inspect_xlsx_schema import exact_source_download,schema_only
CSV_ID=28882431;CSV_BYTES=16409599;CSV_MD5="269d3297531f9d607a0f12a9e3fbf7e4"
def csv_schema():
    url="https://ndownloader.figshare.com/files/"+str(CSV_ID)
    with urllib.request.urlopen(urllib.request.Request(url,
         headers={"User-Agent":"TTF-CaliPopGen-v1-serialization-census"}),timeout=90) as r:
        data=r.read(CSV_BYTES+1)
    if len(data)!=CSV_BYTES or hashlib.md5(data).hexdigest()!=CSV_MD5:
        raise ValueError("Published CSV v1 source identity mismatch")
    with io.TextIOWrapper(io.BytesIO(data),encoding="utf-8-sig",newline="") as fh:
        rdr=csv.reader(fh)
        header=next(rdr)
        rows=0
        widths={}
        for record in rdr:
            rows+=1
            widths[len(record)]=widths.get(len(record),0)+1
    return {"rows":rows,"columns":header,"column_count":len(header),
            "row_width_counts":widths,"md5":CSV_MD5,"GenDist_numeric_values_inspected":False}
if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--source-lock",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    p=a.parse_args()
    lock=json.loads(p.source_lock.read_text())
    x=schema_only(exact_source_download(lock));y=csv_schema()
    out={"schema":"ttf_calipopgen_xlsx_csv_v1_structural_reconciliation",
         "status":"SOURCE_FORMATS_STRUCTURALLY_COMPARED_NO_GENETIC_VALUES",
         "xlsx_rows":x["rows_in_first_worksheet"],"csv_rows":y["rows"],
         "xlsx_columns":x["columns"],"csv_columns":y["columns"],
         "xlsx_column_count":x["column_count"],"csv_column_count":y["column_count"],
         "csv_row_width_counts":y["row_width_counts"],
         "same_number_of_population_pair_records":x["rows_in_first_worksheet"]==y["rows"],
         "Fst_numeric_values_opened":False}
    p.output.parent.mkdir(parents=True,exist_ok=True)
    p.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"xlsx_rows":out["xlsx_rows"],"csv_rows":out["csv_rows"],
       "xlsx_columns":out["xlsx_column_count"],"csv_columns":out["csv_column_count"],
       "same_count":out["same_number_of_population_pair_records"]}))
