"""MultiQC custom module for popmix result aggregation.

This module provides MultiQC integration for population genetics analyses,
including F-statistics, qpAdm, qpWave, and DATES results.

Usage:
    pip install popmix[multiqc]
    multiqc results/ --module popmix
"""

from popmix.multiqc.module import PopmixModule

__all__ = ["PopmixModule"]
