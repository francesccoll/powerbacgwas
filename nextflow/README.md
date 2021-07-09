This directory contains all required files to run PowerBacGWAS as a Nextflow pipeline

## Notes on editing nextflow.config file
* Edit the 'process.container' variable to include the repository/name:tag of the PowerBacGWAS Docker image. By default, the francesccoll/powerbacgwas:amd64 image will be used, which will work on a Linux amd64 host architecture.
