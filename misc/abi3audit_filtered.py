#!/usr/bin/env python3
"""Filtered abi3audit wrapper that allows buffer protocol symbols (stable since 3.0, in Limited API since 3.11)."""
import sys
import json
import subprocess

# Buffer protocol functions: stable since Python 3.0, officially in Limited API since 3.11
ALLOWED_SYMBOLS = {
    'PyObject_GetBuffer',
    'PyBuffer_Release',
    'PyBuffer_FillInfo',
    'PyBuffer_IsContiguous',
}


def main():
    if len(sys.argv) < 2:
        print("Usage: python abi3audit_filtered.py <wheel.whl>")
        return 1

    wheel = sys.argv[1]
    result = subprocess.run(
        ['abi3audit', '--report', wheel],
        capture_output=True, text=True
    )

    try:
        data = json.loads(result.stdout)
    except json.JSONDecodeError:
        print("Failed to parse abi3audit output:")
        print(result.stdout)
        print(result.stderr)
        return 1

    unexpected = []
    allowed_found = []

    for spec, info in data.get('specs', {}).items():
        for so in info.get('wheel', []):
            result_data = so.get('result', {})
            for sym in result_data.get('non_abi3_symbols', []):
                if sym in ALLOWED_SYMBOLS:
                    allowed_found.append(f"{so['name']}: {sym}")
                else:
                    unexpected.append(f"{so['name']}: {sym}")

    if allowed_found:
        print(f"Allowed buffer protocol symbols ({len(allowed_found)}):")
        for v in allowed_found[:10]:  # Show first 10
            print(f"  - {v}")
        if len(allowed_found) > 10:
            print(f"  ... and {len(allowed_found) - 10} more")

    if unexpected:
        print(f"\nERROR: Unexpected abi3 violations ({len(unexpected)}):")
        for v in unexpected:
            print(f"  - {v}")
        return 1

    print("\nabi3audit: OK (only expected buffer protocol symbols)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
