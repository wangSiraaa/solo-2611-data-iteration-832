"""PREPQUANT — exact dilution-ledger and sample-quantitation package.

Two independent components share only rational arithmetic and an exported product
format:

* :class:`prepquant.ledger.Ledger` — ordered, atomic, replayable preparation log
  with provenance books and content-addressed product versions.
* :class:`prepquant.vault.Vault` — measurements pinned to exported products and
  jointly solved exact interval constraints.

All numbers are :class:`fractions.Fraction`; nothing is ever rounded.
"""

from .rational import q, qstr, ZERO
from .expr import LinExpr
from .ledger import Ledger, LedgerError, Product, product_from_payload
from .lp import Infeasible, Unbounded, SolverError
from .vault import Vault, VaultError

__all__ = [
    "q",
    "qstr",
    "ZERO",
    "LinExpr",
    "Ledger",
    "LedgerError",
    "Product",
    "product_from_payload",
    "Vault",
    "VaultError",
    "Infeasible",
    "Unbounded",
    "SolverError",
]
