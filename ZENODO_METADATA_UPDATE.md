# Zenodo v0.2.1 metadata update

Target record: https://zenodo.org/records/21764196

The repository-root `.zenodo.json` is the source of truth for future GitHub-integrated deposits. The existing v0.2.1 record should be edited to match it, particularly:

- replace the current description with the `.zenodo.json` description;
- set the creator affiliation to `University of Michigan`;
- add the eight keywords in `.zenodo.json`;
- add every entry in `.zenodo.json` → `references` to the Zenodo References field;
- retain the MIT license, v0.2.1 version, publication date, release relation, concept DOI, and repository link;
- correct “30 automated checks” to “36 automated checks.”

The references intentionally credit exact engines, scholarly sources, upstream public repositories, and Sagnik Bhadury's development forks. They do not disclose private repositories, working manuscripts, or unpublished method ideas.
