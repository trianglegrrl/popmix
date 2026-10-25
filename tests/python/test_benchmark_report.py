"""The HTML benchmark report (benchmarks/report_html.py, TASKS.md E6)."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2] / "benchmarks"))

import report_html  # noqa: E402

RESULTS = sorted((Path(__file__).parents[2] / "benchmarks" / "results").glob("razorback-*.json"))[-1]


def _payload(page: str) -> dict:
    start = page.index("const DATA = ") + len("const DATA = ")
    return json.loads(page[start:page.index(";\n", start)].replace("<\\/", "</"))


def test_embeds_the_results():
    results = json.loads(RESULTS.read_text())
    page = report_html.render(results)
    data = _payload(page)
    assert data["cases"] == results["cases"]
    assert data["scaling"] is None
    assert page.startswith("<title>")


def test_text_in_the_results_cannot_close_the_script():
    results = {"date": "d", "machine": {}, "versions": {}, "repeats": 1,
               "cases": {"x": {"analysis": "</script><b>"}}}
    page = report_html.render(results, {"curves": {}})
    assert page.count("</script>") == 1
    assert _payload(page)["cases"]["x"]["analysis"] == "</script><b>"
