# Security and data handling

## Public deployment boundary

This application is designed for de-identified research data. A public deployment is not an appropriate destination for protected health information, direct identifiers, secrets, credentials, or data governed by an agreement that prohibits third-party processing.

## Session behavior

- Uploaded files are read into the active R process.
- The application does not intentionally write uploads to a permanent database or repository.
- Temporary files and generated bundles are scoped to the running process/session.
- Session termination is expected to remove temporary files, subject to the host platform's implementation.

Operators must independently configure HTTPS, access control, log retention, container isolation, upload limits, process timeouts, malware scanning where appropriate, and deletion of abandoned temporary storage.

## Reporting vulnerabilities

Report a suspected vulnerability privately to the repository owner. Do not include patient data or confidential datasets in a GitHub issue.
