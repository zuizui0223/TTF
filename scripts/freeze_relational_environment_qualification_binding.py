#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--workflow-run-id",type=int,required=True)
    ap.add_argument("--artifact-id",type=int,required=True)
    ap.add_argument("--head-sha",required=True)
    ap.add_argument("--pregenetic-state",type=Path,required=True)
    ap.add_argument("--geometry-reconstruction",type=Path,required=True)
    ap.add_argument("--opportunity-summary",type=Path,required=True)
    ap.add_argument("--opportunity-design",type=Path,required=True)
    ap.add_argument("--qualification-result",type=Path)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    state=json.loads(args.pregenetic_state.read_text())
    if state.get("schema")!="ttf_relational_environment_pregenetic_state_v0.3":
        raise RuntimeError("unexpected Study-B pregenetic state schema")
    if state.get("response_opened") is not False:
        raise RuntimeError("Study-B pregenetic artifact indicates response opening")
    program_state=state.get("program_state")
    allowed={
        "PASS_TO_CONFIRMATORY_CHARACTER_MASK_PREPARATION",
        "NOT_EVALUABLE_B",
    }
    if program_state not in allowed:
        raise RuntimeError(f"Study-B qualification artifact is incomplete: {program_state!r}")

    files={
        "relational_environment_pregenetic_state_v0.3.json":sha256_path(args.pregenetic_state),
        "relational_environment_geometry_reconstruction_v0.2.json":sha256_path(args.geometry_reconstruction),
        "relational_environment_opportunity_design_v0.2.json":sha256_path(args.opportunity_summary),
        "relational_environment_opportunity_design_v0.2.npz":sha256_path(args.opportunity_design),
    }
    qualification_status=None
    if args.qualification_result is not None and args.qualification_result.is_file():
        q=json.loads(args.qualification_result.read_text())
        if q.get("schema")!="ttf_relational_environment_qualification_result_v0.2":
            raise RuntimeError("unexpected Study-B qualification schema")
        if any(bool(v) for v in q["response_firewall"].values()):
            raise RuntimeError("Study-B qualification response firewall is open")
        qualification_status=q["status"]
        files["relational_environment_qualification_result_v0.2.json"]=sha256_path(args.qualification_result)

    if program_state=="PASS_TO_CONFIRMATORY_CHARACTER_MASK_PREPARATION":
        if qualification_status!="PASS_TO_CONFIRMATORY_CHARACTER_MASK_PREPARATION":
            raise RuntimeError("PASS binding requires a passing qualification result")
    payload={
        "schema":"ttf_relational_environment_qualification_artifact_binding_v0.1",
        "status":"FROZEN_RESPONSE_BLIND_STUDY_B_QUALIFICATION_ARTIFACT",
        "workflow_run_id":int(args.workflow_run_id),
        "artifact_id":int(args.artifact_id),
        "artifact_name":"relational-environment-qualification-v0.2",
        "head_sha":str(args.head_sha),
        "program_state":program_state,
        "qualification_status":qualification_status,
        "files_sha256":files,
        "response_firewall":{
            "Study_B_sequence_identity_opened":False,
            "Study_B_pairwise_genetic_distances_opened":False,
            "Study_B_T_st_computed":False,
            "Study_B_beta_R_computed":False,
        },
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "program_state":program_state,
        "workflow_run_id":payload["workflow_run_id"],
        "artifact_id":payload["artifact_id"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
