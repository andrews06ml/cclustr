## R CMD check results

0 errors | 0 warnings | 2 notes

* checking for hidden files and directories: NOTE
  Found .github directory. This directory contains GitHub Actions 
  workflows for CI/CD and is listed in .Rbuildignore. It will not 
  be included in the package build.

* checking for future file timestamps: NOTE
  unable to verify current time. This is a server-side issue on the 
  check platform, not related to the package.

## Test environments

### Local
* Windows 11, R 4.3.3: 0 errors | 0 warnings | 2 notes

### rhub
* Ubuntu 24.04 LTS, R-devel (2026-04-26 r89963): 0 errors | 0 warnings | 1 note
* Windows Server 2022, R-devel (2026-04-26 r89963): 0 errors | 0 warnings | 1 note
* Fedora Linux 42 / ATLAS, R-devel (2026-04-26 r89963): 0 errors | 0 warnings | 1 note

### win-builder
* Windows, R-devel (2026-04-27 r89967): 0 errors | 0 warnings | 2 notes
* Windows, R-release (R 4.6.0): 0 errors | 0 warnings | 2 notes

## Notes

* This is a new submission.
* The .github directory is excluded from the package build via
  .Rbuildignore.
* The "unable to verify current time" note is a server-side issue
  on the check platform, not related to the package.
