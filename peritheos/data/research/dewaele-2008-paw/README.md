# Dewaele 2008 PAW research collection

54 published Table VI PAW parameterizations in nine native EOSMAT files.
Research only: excluded from the main material database and catalog APIs.
The 13 experimental Table IV EOS remain in the main catalog.

All original coefficients, errors, source and validation findings are retained.
Figure 4 graphical checkpoints are in datasets/; complete checks are in
reproduction.json, primary-source findings in source-audit.json, and the
publication transcription in publication-source.json. Fifty sets have
source-equation parity; two bcc-Fe LDA sets are deferred and two hcp-Fe
LDA sets retain their unresolved semicore-identity finding.

Open an EOSMAT file explicitly with load_eosmat and Material.from_eosmat.
Load its dataset metadata with Dataset.from_mapping(metadata,
resource_root=collection_directory). The companion documentation is
docs/research/dewaele-2008-paw.md in the repository.

Primary source: https://doi.org/10.1103/PhysRevB.78.104102
Source numerical facts retain attribution; the publisher PDF is not
redistributed and no licence over third-party content is asserted.
