# Windows release signature audit

For [issue #4305](https://github.com/casadi/casadi/issues/4305), run the
**Windows signature audit** workflow with a stable release tag (default: 3.8.0).
It downloads the published CPython 3.11+ x64 wheel and MATLAB x64 ZIP from
`casadi/casadi`, records their SHA256 hashes, inventories Authenticode signatures
on native binaries, and tests the installed release rather than a source build.

Python tests run in separate processes for symbolic evaluation/AD, IPOPT,
CVODES, JIT, and a compiled external function. MATLAB tests exercise the MEX,
symbolic evaluation/AD, and IPOPT. Generated DLLs and per-case Python logs are
retained; MATLAB output is in the Actions log.

The runner is Windows Server 2025, **not consumer Windows 11 with Smart App
Control enabled**. A passing workflow is not SAC certification and does not
answer whether Microsoft currently trusts a particular unsigned file's hash.

After baseline tests, the workflow attempts to install Microsoft's published
`SmartAppControlAuditNoISG.bin` as an additional audit-only App Control policy.
It does not turn on enforcement or change the runner's SAC registry settings.
The policy archive hash is pinned. Installation failure is reported as unavailable
coverage; baseline inventory and tests remain useful. Even successful installation
does not establish coverage: the report requires a Code Integrity event 3076
from that policy for the locally compiled `unsigned_probe.dll` control.

Read `summary.md` and `coverage.json` in the uploaded artifact first. Full
signature JSON/CSV, runner/policy metadata, Code Integrity events (including raw
XML), artifact hashes, and test outputs provide supporting evidence. Unsigned
files are findings, not workflow failures. Smoke-test failures fail the workflow.

The first implementation step after reviewing the inventory is trusted RSA
Authenticode signing of final packaged binaries and dependencies before ZIP/wheel
creation, followed by signature verification. Existing vendor signatures should
be preserved. Locally generated JIT/external DLLs require separate investigation.

References:

- [Microsoft SAC testing and audit policies](https://learn.microsoft.com/en-us/windows/apps/develop/smart-app-control/test-your-app-with-smart-app-control)
- [Microsoft SAC signing requirements](https://learn.microsoft.com/en-us/windows/apps/develop/smart-app-control/code-signing-for-smart-app-control)
- [GitHub Windows runner image](https://github.com/actions/runner-images/blob/main/images/windows/Windows2025-Readme.md)
