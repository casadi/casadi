param(
    [Parameter(Mandatory)]
    [ValidateSet('Inventory', 'InstallPolicy', 'Collect')]
    [string] $Action
)

$ErrorActionPreference = 'Stop'
$policyId = '5283AC0F-FFF1-49AE-ADA1-8A933130CAD6'
$reports = Join-Path $PWD 'reports'
New-Item -ItemType Directory -Force $reports | Out-Null

function Export-Signatures([string] $Root, [string] $Name) {
    $rootPath = (Resolve-Path $Root).Path
    $rows = @(foreach ($file in Get-ChildItem $rootPath -Recurse -File) {
        # Inspect the PE header rather than relying on extensions (MEX, PYD, DLL).
        $stream = $file.OpenRead()
        try { $isPE = $stream.ReadByte() -eq 0x4d -and $stream.ReadByte() -eq 0x5a }
        finally { $stream.Dispose() }
        if (-not $isPE) { continue }
        $signature = Get-AuthenticodeSignature -LiteralPath $file.FullName
        [pscustomobject]@{
            Package = $Name
            Path = [IO.Path]::GetRelativePath($rootPath, $file.FullName)
            SHA256 = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash
            Status = [string]$signature.Status
            StatusMessage = $signature.StatusMessage
            SignatureType = [string]$signature.SignatureType
            Subject = $signature.SignerCertificate.Subject
            Thumbprint = $signature.SignerCertificate.Thumbprint
            PublicKeyAlgorithm = $signature.SignerCertificate.PublicKey.Oid.Value
            TimestampSubject = $signature.TimeStamperCertificate.Subject
        }
    })
    if ($rows.Count -eq 0) { throw "No PE binaries found under $Root" }
    $rows | ConvertTo-Json -Depth 5 -AsArray | Set-Content "$reports/$Name-signatures.json"
    $rows | Export-Csv "$reports/$Name-signatures.csv" -NoTypeInformation
    $rows | Group-Object Status | Select-Object Name, Count | Format-Table
}

switch ($Action) {
    'Inventory' {
        $os = Get-CimInstance Win32_OperatingSystem
        $sac = Get-ItemProperty 'HKLM:\SYSTEM\CurrentControlSet\Control\CI\Policy' -ErrorAction SilentlyContinue
        [pscustomobject]@{
            OS = $os.Caption
            Version = $os.Version
            Build = $os.BuildNumber
            ImageOS = $env:ImageOS
            ImageVersion = $env:ImageVersion
            Release = $env:CASADI_RELEASE
            RunURL = "https://github.com/$env:GITHUB_REPOSITORY/actions/runs/$env:GITHUB_RUN_ID"
            Commit = $env:GITHUB_SHA
            VerifiedAndReputablePolicyState = $sac.VerifiedAndReputablePolicyState
            Limitation = 'Windows Server signature audit and smoke tests; not Windows 11 consumer SAC enforcement or cloud-reputation certification.'
        } | ConvertTo-Json | Set-Content "$reports/environment.json"
        if (Get-Command CiTool.exe -ErrorAction SilentlyContinue) {
            & CiTool.exe -lp -json | Set-Content "$reports/policies-before.json"
            if ($LASTEXITCODE) { throw 'Could not enumerate Code Integrity policies' }
        }
        Export-Signatures $env:CASADI_PACKAGE 'python'
        Export-Signatures 'matlab-package' 'matlab'
    }
    'InstallPolicy' {
        # Published by Microsoft, linked from the SAC developer testing guide.
        # Audit only: no enforcement, no changes to existing policies, no registry SAC override.
        $url = 'https://download.microsoft.com/download/b/4/5/b45e7463-6ae0-461d-95ff-89cec7ce5159/SAC%20Audit%20Policies.zip'
        Invoke-WebRequest $url -OutFile "$reports/policies.zip"
        $hash = (Get-FileHash "$reports/policies.zip" -Algorithm SHA256).Hash
        if ($hash -ne 'D7D729412DA6D1B759B460043B401961EF125F1B194574DFFF1F53DD985B07C7') {
            throw 'Microsoft audit policy archive changed; review before updating its pinned hash'
        }
        Expand-Archive "$reports/policies.zip" -DestinationPath "$reports/policies" -Force
        $policy = "$reports/policies/{$policyId}.cip"
        Copy-Item "$reports/policies/SmartAppControlAuditNoISG.bin" $policy
        (Get-Date).ToUniversalTime().ToString('o') | Set-Content "$reports/audit-start.txt"
        & wevtutil.exe sl Microsoft-Windows-CodeIntegrity/Operational /e:true
        if ($LASTEXITCODE) { throw 'Could not enable Code Integrity event collection' }
        & CiTool.exe -up $policy -json | Tee-Object "$reports/policy-install.json"
        if ($LASTEXITCODE) { throw 'Runner cannot install Microsoft signing-only audit policy' }
        $policies = & CiTool.exe -lp -json
        $policies | Set-Content "$reports/policies-after.json"
        if ($LASTEXITCODE -or ($policies -join "`n") -notmatch $policyId) {
            throw 'The requested audit policy was not found in the installed policy list'
        }
        # A successful update alone does not prove auditing works. Collect requires
        # an event from this policy for our newly compiled unsigned probe DLL.
    }
    'Collect' {
        $events = @()
        $eventError = $null
        if (Test-Path "$reports/audit-start.txt") {
            $start = [datetime]::Parse((Get-Content "$reports/audit-start.txt" -Raw))
            # Flush the asynchronous event pipeline before collecting evidence.
            Start-Sleep -Seconds 5
            try {
                $events = @(Get-WinEvent -FilterHashtable @{
                    LogName = 'Microsoft-Windows-CodeIntegrity/Operational'
                    StartTime = $start
                    Id = 3076, 3077, 3089, 3099
                } -ErrorAction Stop | ForEach-Object {
                    [xml]$xml = $_.ToXml()
                    $data = @{}
                    foreach ($field in $xml.Event.EventData.Data) {
                        $data[[string]$field.Name] = [string]$field.'#text'
                    }
                    [pscustomobject]@{
                        Id = $_.Id
                        TimeCreated = $_.TimeCreated.ToUniversalTime().ToString('o')
                        Data = $data
                        Message = $_.Message
                        XML = $_.ToXml()
                    }
                })
            } catch {
                $eventError = $_.Exception.Message
                Write-Warning "Code Integrity event query: $eventError"
            }
        }
        ConvertTo-Json -InputObject $events -Depth 8 | Set-Content "$reports/code-integrity-events.json"
        $policyEvents = @($events | Where-Object {
            $_.Id -eq 3076 -and $_.XML -match $policyId
        })
        $probeEvents = @($policyEvents | Where-Object { $_.XML -match 'unsigned_probe\.dll' })
        $casadiEvents = @($policyEvents | Where-Object {
            $_.XML -match '(?i)(site-packages\\casadi\\|matlab-package\\|reports\\audit\\)'
        })
        ConvertTo-Json -InputObject $casadiEvents -Depth 8 | Set-Content "$reports/casadi-audit-events.json"
        if (Test-Path "$reports/baseline/external/unsigned_probe.dll") {
            Export-Signatures "$reports/baseline" 'generated-baseline'
        }
        if (Test-Path "$reports/audit/external/unsigned_probe.dll") {
            Export-Signatures "$reports/audit" 'generated-audit'
        }
        $coverage = [pscustomobject]@{
            PolicyStep = $env:POLICY_OUTCOME
            AuditProbeObserved = $probeEvents.Count -gt 0
            PolicyAuditEvents = $policyEvents.Count
            CasadiAuditEvents = $casadiEvents.Count
            EventQueryError = $eventError
            PythonBaseline = $env:PYTHON_BASELINE
            MatlabBaseline = $env:MATLAB_BASELINE
            PythonAudit = $env:PYTHON_AUDIT
            MatlabAudit = $env:MATLAB_AUDIT
            ConsumerSACEnforcementTested = $false
            CloudReputationTested = $false
        }
        $coverage | ConvertTo-Json | Set-Content "$reports/coverage.json"
        $summary = @(
            "## Windows signature audit: CasADi $env:CASADI_RELEASE"
            ''
            '**Scope:** Windows Server signature inventory and runtime smoke tests. Consumer Windows 11 SAC enforcement and cloud reputation were not tested.'
            ''
            '| Package | Signature status | Native files |'
            '| --- | --- | ---: |'
        )
        foreach ($name in 'python', 'matlab', 'generated-baseline', 'generated-audit') {
            if (Test-Path "$reports/$name-signatures.json") {
                $rows = Get-Content "$reports/$name-signatures.json" -Raw | ConvertFrom-Json
                foreach ($group in $rows | Group-Object Status) {
                    $summary += "| $name | $($group.Name) | $($group.Count) |"
                }
            }
        }
        $summary += @(
            ''
            "Python baseline: **$env:PYTHON_BASELINE**. MATLAB baseline: **$env:MATLAB_BASELINE**."
            ''
            "Audit policy installation: **$env:POLICY_OUTCOME**. Unsigned probe audit event observed: **$($coverage.AuditProbeObserved)**."
            ''
            "Python under requested audit policy: **$env:PYTHON_AUDIT**. MATLAB: **$env:MATLAB_AUDIT**."
            ''
            "CasADi/generated-code events from the signing-only policy: **$($casadiEvents.Count)**."
            ''
            'Unsigned files lack the signature fallback when reputation is unknown. This does not establish that consumer SAC currently blocks those files. Audit-mode execution success does not mean the file would be allowed in enforcement mode.'
        )
        if (-not $coverage.AuditProbeObserved) {
            $summary += @('', '**Signing-only policy coverage was NOT demonstrated.** Use the signature inventory and baseline results only; inspect policy/event diagnostics in the artifact.')
        }
        $summary | Set-Content "$reports/summary.md"
        $summary | Write-Output
        if ($env:GITHUB_STEP_SUMMARY) { $summary >> $env:GITHUB_STEP_SUMMARY }
    }
}
