"""Main MultiQC module class for popmix.

This module provides the entry point for MultiQC integration.
"""

from typing import Any


class PopmixModule:
    """MultiQC module for popmix population genetics results.

    This module aggregates and displays:
    - F-statistics (f2, f3, f4, FST)
    - qpAdm admixture weights
    - qpWave rank tests
    - DATES admixture timing

    Attributes
    ----------
    name : str
        Module name for MultiQC
    anchor : str
        HTML anchor for navigation
    info : str
        Brief description
    """

    name = "popmix"
    anchor = "popmix"
    info = "Population genetics F-statistics and admixture analysis"

    def __init__(self):
        """Initialize the popmix MultiQC module."""
        self.data: dict[str, dict[str, Any]] = {}
        self.qpadm_data: list[dict[str, Any]] = []
        self.dates_data: list[dict[str, Any]] = []
        self.fstat_data: list[dict[str, Any]] = []

    @classmethod
    def get_name(cls) -> str:
        """Get module name."""
        return cls.name

    @classmethod
    def get_anchor(cls) -> str:
        """Get HTML anchor."""
        return cls.anchor

    def find_log_files(self, directory: str) -> list[str]:
        """Find popmix output files in a directory.

        Parameters
        ----------
        directory : str
            Directory to search

        Returns
        -------
        List[str]
            List of found file paths
        """
        from popmix.multiqc.patterns import find_popmix_files
        return find_popmix_files(directory)

    def parse_logs(self, files: list[str]) -> None:
        """Parse found log files.

        Parameters
        ----------
        files : List[str]
            List of file paths to parse
        """
        from popmix.multiqc.parsers import parse_popmix_json
        from popmix.multiqc.patterns import get_sample_name

        for filepath in files:
            try:
                data = parse_popmix_json(filepath)
                sample_name = get_sample_name(filepath)

                result_type = data.get("type", "unknown")

                if result_type == "qpadm":
                    self.qpadm_data.append({"sample": sample_name, **data})
                elif result_type == "dates":
                    self.dates_data.append({"sample": sample_name, **data})
                elif result_type in ("f2", "f3", "f4", "fst"):
                    self.fstat_data.append({"sample": sample_name, **data})

                # Store for general stats
                self.data[sample_name] = data

            except Exception as e:
                # Log error but continue processing
                print(f"Warning: Could not parse {filepath}: {e}")

    def get_general_stats_data(self) -> dict[str, dict[str, Any]]:
        """Get data for General Statistics table.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            Data formatted for MultiQC general stats
        """
        from popmix.multiqc.general_stats import format_general_stats_data
        return format_general_stats_data(self.data)

    def get_general_stats_headers(self) -> dict[str, dict[str, Any]]:
        """Get General Statistics column headers.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            Column header definitions
        """
        from popmix.multiqc.general_stats import get_general_stats_headers
        return get_general_stats_headers()

    def generate_sections(self) -> list[dict[str, Any]]:
        """Generate report sections.

        Returns
        -------
        List[Dict[str, Any]]
            List of section definitions
        """
        sections = []

        if self.qpadm_data:
            sections.append({
                "name": "qpAdm Results",
                "anchor": "popmix-qpadm",
                "content": self._generate_qpadm_section(),
            })

        if self.dates_data:
            sections.append({
                "name": "DATES Results",
                "anchor": "popmix-dates",
                "content": self._generate_dates_section(),
            })

        if self.fstat_data:
            sections.append({
                "name": "F-Statistics",
                "anchor": "popmix-fstats",
                "content": self._generate_fstat_section(),
            })

        return sections

    def _generate_qpadm_section(self) -> str:
        """Generate qpAdm section HTML."""
        if not self.qpadm_data:
            return "<p>No qpAdm results found.</p>"

        rows = []
        for data in self.qpadm_data:
            target = data.get("target", "")
            sources = ", ".join(data.get("sources", []))
            p_value = data.get("p_value", 0)
            fit_class = "good" if p_value > 0.05 else "poor"

            rows.append(f"""
            <tr>
                <td>{target}</td>
                <td>{sources}</td>
                <td>{p_value:.4f}</td>
                <td class="{fit_class}">{fit_class.title()}</td>
            </tr>
            """)

        return f"""
        <table class="table table-bordered">
            <thead>
                <tr>
                    <th>Target</th>
                    <th>Sources</th>
                    <th>P-value</th>
                    <th>Fit</th>
                </tr>
            </thead>
            <tbody>
                {"".join(rows)}
            </tbody>
        </table>
        """

    def _generate_dates_section(self) -> str:
        """Generate DATES section HTML."""
        if not self.dates_data:
            return "<p>No DATES results found.</p>"

        rows = []
        for data in self.dates_data:
            pop = data.get("test_pop", "")
            gens = data.get("generations", 0)
            years = data.get("years_bp", 0)
            r2 = data.get("r_squared", 0)

            rows.append(f"""
            <tr>
                <td>{pop}</td>
                <td>{gens:.0f}</td>
                <td>{years:,.0f}</td>
                <td>{r2:.3f}</td>
            </tr>
            """)

        return f"""
        <table class="table table-bordered">
            <thead>
                <tr>
                    <th>Population</th>
                    <th>Generations</th>
                    <th>Years BP</th>
                    <th>R²</th>
                </tr>
            </thead>
            <tbody>
                {"".join(rows)}
            </tbody>
        </table>
        """

    def _generate_fstat_section(self) -> str:
        """Generate F-statistics section HTML."""
        if not self.fstat_data:
            return "<p>No F-statistic results found.</p>"

        rows = []
        for data in self.fstat_data:
            stat_type = data.get("type", "").upper()
            estimate = data.get("estimate", 0)
            z = data.get("z_score", 0)
            sig = "Yes" if abs(z) >= 2 else "No"

            rows.append(f"""
            <tr>
                <td>{stat_type}</td>
                <td>{estimate:.6f}</td>
                <td>{z:.2f}</td>
                <td>{sig}</td>
            </tr>
            """)

        return f"""
        <table class="table table-bordered">
            <thead>
                <tr>
                    <th>Statistic</th>
                    <th>Estimate</th>
                    <th>Z-score</th>
                    <th>Significant</th>
                </tr>
            </thead>
            <tbody>
                {"".join(rows)}
            </tbody>
        </table>
        """
