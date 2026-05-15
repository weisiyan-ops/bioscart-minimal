"""Minimal BioSCART research prototype.

This package implements the first CT-only workflow:

- import a DICOM CT series and RTSTRUCT,
- find a GTV ROI,
- rasterize GTV to the CT grid,
- create five transparent research subregions,
- optionally extract PyRadiomics features,
- emit research-only recommendations.
"""

__version__ = "0.1.0"

