# Powershell script for running multiple FACTS files present in a given directory using command-line FACTS
$Dest = '..\RegressionTests\Tests\DF\TimeToEvent F*' # Specify the directory containing the FACTS files to run here.

Write-Host "Started FACTS tests run"

Get-ChildItem -Path $Dest -Recurse |
Where-Object {$_.Name -match '_results'} |
Remove-Item -Recurse

$FactsFileFullNames = Get-ChildItem -Include '*.facts' -Path $Dest -Recurse | select -ExpandProperty FullName

$Process = "C:\Program Files\BerryConsultants\FACTS 6.4.0\FACTS File Loader.exe"

foreach ($i in $FactsFileFullNames) {
Start-Process $Process -ArgumentList "-b -f `"$i`" -n 4 -p 2 -nSubjectFiles 4" -Wait -NoNewWindow
}

Write-Host "Completed FACTS tests run"