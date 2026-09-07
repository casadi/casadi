# CasADi 3.8.0 Windows audit findings

The [Windows CI audit](https://github.com/casadi/casadi/actions/runs/34151295106)
completed successfully on 2026-09-07, testing commit
`0cdc4d76b083e32503cacbdaf1282279fd66530e`.

The runner was Windows Server 2025 Datacenter, build 26100, image
`win25-vs2026` / `20260824.214.3`, with Python 3.11.9 and MATLAB R2024b.
The tested binaries came from the public 3.8.0 release, not this branch.

## Signature inventory

| Package | Native binaries | Windows Authenticode result |
| --- | ---: | --- |
| CPython 3.11+ x64 wheel | 117 | All `NotSigned` |
| MATLAB x64 ZIP | 144 | All `NotSigned` |
| Locally generated JIT and external DLLs, per test phase | 2 | Both `NotSigned` |

An independent inspection of the two release archives also found no embedded
PE certificate tables. Windows' signature inventory agrees with that result.

Artifact SHA256 hashes:

- `casadi-3.8.0-cp311-abi3-win_amd64.whl`:
  `8630f9c68c7d05b49b538f05cc4f0873f0eaac720153c324f3dd00fc222b4d79`
- `casadi-3.8.0-windows64-matlab2018b.zip`:
  `e1b9af1a0557e26fffbfce341c20e9d5482ce8351539f309559e32c2a7acc20e`

## Runtime and policy experiment

All five Python cases passed both before and after installing the audit policy:

- Symbolic evaluation and automatic differentiation.
- An IPOPT solve with a checked solution.
- CVODES integration checked against an analytic solution.
- JIT compilation and evaluation.
- Compiling generated C with MSVC and loading/evaluating the external DLL.

MATLAB MEX loading, symbolic evaluation/AD, and IPOPT also passed in both phases.

Microsoft's `VerifiedAndReputableDesktopEvaluationAuditNoISG` policy installed
and became active with `Enabled:Audit Mode`. Its policy GUID was
`{5283ac0f-fff1-49ae-ada1-8a933130cad6}`. CiTool's `IsEnforced: true` means
that this policy is active; its explicit audit-mode option means it logs
would-block events instead of blocking execution.

The policy produced **42 events with ID 3076**, including **33 events for
CasADi or generated code**, covering 31 distinct paths. These included:

- `_casadi.pyd` and `casadiMEX.mexw64`.
- `libcasadi.dll`, the IPOPT plugin, and `libipopt-3.dll`.
- Bundled OpenBLAS, MUMPS, METIS, GCC, Fortran, and threading runtime DLLs.
- The CVODES and shell compiler plugins.
- Both the JIT DLL and `unsigned_probe.dll`, our unsigned control.

The control event verifies that the signing-only audit actually ran. Passing
runtime tests under an audit policy does not mean the files would be allowed
under the corresponding enforcement policy.

## Interpretation and next change

This establishes a release-signing gap across the bindings and their native
runtime dependencies. It does **not** answer whether consumer Smart App Control
currently trusts these hashes through its cloud reputation service. Windows 11
consumer SAC enforcement and cloud reputation were not tested.

The next release-engineering change should sign the final Windows native payload
with a publicly trusted RSA Authenticode certificate and timestamp, then verify
all expected signatures before publication. In the existing binary workflow,
Python signing belongs after the last dependency DLL copy and before ZIP/wheel
creation; MATLAB needs the equivalent final-payload signing step. Preserve valid
vendor signatures when present. Wheel RECORD hashes must be generated after
signing. Re-run this audit against the resulting release candidate.

JIT and user-compiled external DLLs are a separate compatibility problem: signing
our distributed files cannot sign binaries generated later on the user's machine.

Full JSON/CSV inventories, raw Code Integrity events, policy state, generated
DLLs, and Python logs are in the run's `windows-signature-audit-3.8.0` artifact.
MATLAB output is in the job log.
