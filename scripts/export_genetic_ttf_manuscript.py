#!/usr/bin/env python3
"""Export validated scalar receipts without opening any genetic source files."""
from __future__ import annotations

import argparse
from pathlib import Path

from ttf.genetic_manuscript_export import build_closed_export, build_empirical_export, write_export


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--result', type=Path, help='Completed Phase-4 result with authorization hash')
    mode.add_argument('--closed-receipt', type=Path, help='Canonical v0.3 opening state or failed Phase-1/2/3 receipt')
    parser.add_argument('--authorization', type=Path)
    parser.add_argument('--qualification', type=Path)
    parser.add_argument('--self-qualification', type=Path)
    parser.add_argument('--phase4-rule', type=Path, default=Path('docs/supporting/genetic_phylogatr_phase4_response_rule_v0.1.json'))
    parser.add_argument('--self-rule', type=Path, default=Path('docs/supporting/genetic_phylogatr_phase3_self_detectability_rule_v0.1.json'))
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    supporting = (args.authorization, args.qualification, args.self_qualification)
    if args.result and not all(supporting):
        parser.error('--result requires --authorization, --qualification and --self-qualification')
    if args.closed_receipt and any(supporting):
        parser.error('--closed-receipt cannot be combined with empirical supporting receipts')
    try:
        report = (build_closed_export(args.closed_receipt) if args.closed_receipt else
                  build_empirical_export(result=args.result, authorization=args.authorization,
                      qualification=args.qualification, self_qualification=args.self_qualification,
                      phase4_rule=args.phase4_rule, self_rule=args.self_rule))
        output = write_export(report, args.output_dir)
    except (ValueError, KeyError, TypeError, OSError) as exc:
        parser.error(str(exc))
    print(f"{report.manifest['mode']}: {output}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
