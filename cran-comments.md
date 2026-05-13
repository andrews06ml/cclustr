## Test results
0 errors | 0 warnings | 1 note

- Local Windows 11 installation, R 4.5.2
- win-builder (devel and release)
- R-hub Ubuntu Linux 24.04.4 LTS (R-devel)
- R-hub Fedora Linux 42 (Atlas) (R-devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

* Checking spelling in DESCRIPTION... NOTE
  Possibly misspelled words in DESCRIPTION:
    Monti, Pihur, Senbabaoglu (Author names and specific references)
    et, al (Standard Latin abbreviations for citations)
  
  Note: These terms are verified as author names or correct technical 
  abbreviations. They have been added to the WORDLIST to ensure consistency.

## Resubmission
This is the first submission of the 'cclustr' package.

## Additional Information
- The package maintains compatibility across Windows and Linux (Ubuntu and Fedora).
- Documentation URLs were checked with urlchecker::url_check() and are all functional.
- All examples and vignettes were tested, including those marked as '--run-donttest'.
- The .github directory is excluded from the package build via
  .Rbuildignore.
