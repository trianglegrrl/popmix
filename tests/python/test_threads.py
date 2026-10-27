"""popmix uses the threads it is told to, 4 by default (TASKS.md F7)."""

import subprocess
import sys


def _run(code: str) -> str:
    return subprocess.run([sys.executable, "-c", "import popmix as pm\n" + code],
                          check=True, capture_output=True, text=True).stdout.strip()


def test_the_default_is_four_threads_and_asking_does_not_fix_it():
    assert _run("print(pm.get_threads())\nprint(pm.init_threads(8))\nprint(pm.get_threads())").split() == \
        ["4", "True", "8"]


def test_the_thread_count_can_change_after_a_computation():
    code = ("import os\n"
            "p = os.path.join(os.path.dirname(pm.__file__), '..', '..', 'tests', 'fixtures', 'small_eigenstrat')\n"
            "pm.init_threads(2)\n"
            "pm.extract_f2(p, [x for x, _ in pm.list_populations(p)][:3])\n"
            "print(pm.init_threads(3), pm.get_threads())")
    assert _run(code) == "True 3"


def test_long_computations_release_the_gil():
    """Python threads keep running while popmix computes (TASKS.md F6)."""
    import threading
    import time
    from pathlib import Path

    import popmix as pm

    prefix = str(Path(__file__).parents[1] / "golden" / "dates" / "family_chr20_22")
    pops = [p for p, _ in pm.list_populations(prefix)]
    took = []

    def work() -> None:
        start = time.perf_counter()
        pm.smartpca(prefix, pops, [])
        took.append(time.perf_counter() - start)

    worker = threading.Thread(target=work)
    worker.start()
    # with the GIL held throughout, this loop would run once or twice; a
    # loaded machine delays some of its turns, but not most of them
    ticks = 0
    while worker.is_alive():
        time.sleep(0.001)
        ticks += 1
    worker.join()
    assert took[0] > 0.05, "the computation is too short to tell"
    assert ticks >= 10
