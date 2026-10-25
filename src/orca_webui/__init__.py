"""orca-webui: a Gradio web UI for simple ORCA computational-chemistry calculations.

Public entry point is :func:`orca_webui.app.main` (exposed as the ``orca-webui``
console command). Submodules are imported lazily by that entry point so that
lightweight consumers (e.g. importing :mod:`orca_webui.utils` for its parsing/spectrum
helpers) do not pull in the full Gradio/FastAPI UI stack.
"""

__version__ = "0.1.12"
